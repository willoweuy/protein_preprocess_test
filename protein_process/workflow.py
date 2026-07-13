from __future__ import annotations

from pathlib import Path
import tempfile

from .models import MissingGap, MutationRequest, ProcessingOptions, ProcessingResult
from .pdb_text import (
    PROTONATION_RESNAMES,
    STANDARD_AA,
    add_terminal_caps,
    apply_manual_residue_names,
    apply_split_segids,
    apply_zinc_coordination_names,
    insert_gromacs_ter_records,
    pdb_to_cif_text,
    remove_first_nucleic_phosphate_atoms,
    restore_chain_ids_from_mmcif,
)
from .propka import run_propka_prediction

AA_1_TO_3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}


class LongGapDecisionRequired(RuntimeError):
    def __init__(self, gaps: tuple[MissingGap, ...]):
        super().__init__("Missing internal residue gap exceeds the configured threshold.")
        self.gaps = gaps


def _import_openmm_stack():
    try:
        from openmm.app import PDBFile, PDBxFile
        from pdbfixer import PDBFixer
        from pdbfixer.pdbfixer import Sequence
    except ImportError as exc:
        raise RuntimeError(
            "This module needs `pdbfixer` and `openmm`. Create the conda environment from "
            "`environment.yml`, then run `streamlit run app.py` inside that environment."
        ) from exc
    return PDBFixer, PDBFile, PDBxFile, Sequence


def _format_residue_name(residue) -> str:
    return getattr(residue, "name", str(residue)).upper()


def _chain_id_from_index(fixer, chain_index: int) -> str:
    chains = list(fixer.topology.chains())
    if chain_index >= len(chains):
        return ""
    return chains[chain_index].id or ""


def _residue_id_before_gap(fixer, chain_index: int, residue_index: int) -> str:
    chains = list(fixer.topology.chains())
    if chain_index >= len(chains):
        return ""
    residues = list(chains[chain_index].residues())
    if not residues:
        return ""
    before_index = max(0, min(residue_index - 1, len(residues) - 1))
    return getattr(residues[before_index], "id", "") or ""


def _missing_gaps_from_fixer(fixer) -> tuple[MissingGap, ...]:
    gaps: list[MissingGap] = []
    for (chain_index, residue_index), residues in getattr(fixer, "missingResidues", {}).items():
        names = tuple(_format_residue_name(residue) for residue in residues)
        gaps.append(
            MissingGap(
                chain_id=_chain_id_from_index(fixer, chain_index),
                after_residue=_residue_id_before_gap(fixer, chain_index, residue_index),
                length=len(residues),
                residue_names=names,
            )
        )
    return tuple(gaps)


def _long_gaps(gaps: tuple[MissingGap, ...], threshold: int) -> tuple[MissingGap, ...]:
    return tuple(gap for gap in gaps if gap.length > threshold)


def _remove_long_missing_residue_entries(fixer, threshold: int) -> None:
    fixer.missingResidues = {
        key: residues
        for key, residues in fixer.missingResidues.items()
        if len(residues) <= threshold
    }


def _remove_terminal_missing_residue_entries(fixer) -> list[str]:
    chains = list(fixer.topology.chains())
    kept = {}
    report: list[str] = []
    for (chain_index, residue_index), residues in getattr(fixer, "missingResidues", {}).items():
        chain_id = _chain_id_from_index(fixer, chain_index)
        chain_residues = list(chains[chain_index].residues()) if chain_index < len(chains) else []
        is_terminal = residue_index == 0 or residue_index >= len(chain_residues)
        if is_terminal:
            terminus = "N-terminal" if residue_index == 0 else "C-terminal"
            names = ", ".join(_format_residue_name(residue) for residue in residues) or "unknown"
            report.append(
                f"Ignored {terminus} missing residues for chain {chain_id or '_'}: {len(residues)} residue(s) {names}."
            )
        else:
            kept[(chain_index, residue_index)] = residues
    fixer.missingResidues = kept
    return report


def _topology_residue_lookup(fixer) -> dict[tuple[str, str], str]:
    lookup: dict[tuple[str, str], str] = {}
    for chain in fixer.topology.chains():
        for residue in chain.residues():
            lookup[(chain.id or "", residue.id)] = residue.name.upper()
    return lookup


def _apply_pdbfixer_mutations(
    fixer,
    mutations: tuple[MutationRequest, ...],
) -> tuple[tuple[MutationRequest, ...], list[str]]:
    deferred: list[MutationRequest] = []
    report: list[str] = []
    lookup = _topology_residue_lookup(fixer)
    by_chain: dict[str, list[str]] = {}

    for mutation in mutations:
        new_resname = mutation.new_resname.strip().upper()
        if not new_resname:
            continue
        if new_resname in PROTONATION_RESNAMES or new_resname not in STANDARD_AA:
            deferred.append(mutation)
            continue
        key = (mutation.chain_id.strip(), str(mutation.residue_number).strip())
        old_resname = lookup.get(key)
        if not old_resname:
            deferred.append(mutation)
            report.append(
                f"PDBFixer mutation deferred: chain {key[0]} residue {key[1]} was not found before repair."
            )
            continue
        by_chain.setdefault(key[0], []).append(f"{old_resname}-{key[1]}-{new_resname}")

    for chain_id, mutation_specs in by_chain.items():
        try:
            fixer.applyMutations(mutation_specs, chain_id)
            report.append(f"PDBFixer applied mutation(s) on chain {chain_id}: {', '.join(mutation_specs)}.")
        except Exception as exc:
            report.append(f"PDBFixer could not apply mutation(s) on chain {chain_id}: {exc}")
            for spec in mutation_specs:
                _old, residue_number, new_resname = spec.split("-")
                deferred.append(
                    MutationRequest(chain_id=chain_id, residue_number=residue_number, new_resname=new_resname)
                )

    return tuple(deferred), report


def _write_fixer_to_pdb_text(fixer, PDBFile) -> str:
    with tempfile.NamedTemporaryFile("w+", suffix=".pdb", delete=True) as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
        handle.flush()
        handle.seek(0)
        return handle.read()


def _load_fasta_sequences(fasta_path: Path, fixer, Sequence) -> tuple[list, list[str]]:
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise RuntimeError("Biopython is required for FASTA sequence input.") from exc

    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if not records:
        return [], [f"FASTA sequence file {fasta_path} contained no records."]

    chains = [chain for chain in fixer.topology.chains() if len(list(chain.residues())) > 0]
    chain_ids = [chain.id or "" for chain in chains]
    sequences = []
    report: list[str] = []

    for index, record in enumerate(records):
        record_id = str(record.id).strip()
        chain_id = record_id
        if chain_id not in chain_ids:
            if len(records) == 1 and len(chain_ids) == 1:
                chain_id = chain_ids[0]
                report.append(
                    f"FASTA record {record_id!r} assigned to only chain {chain_id or '_'}."
                )
            elif index < len(chain_ids) and record_id not in chain_ids:
                chain_id = chain_ids[index]
                report.append(
                    f"FASTA record {record_id!r} assigned by order to chain {chain_id or '_'}."
                )
            else:
                report.append(
                    f"FASTA record {record_id!r} skipped because no matching chain ID was found."
                )
                continue

        residues: list[str] = []
        unknown: list[str] = []
        for aa in str(record.seq).upper():
            if aa in {"-", " ", "\n", "\r", "\t"}:
                continue
            resname = AA_1_TO_3.get(aa)
            if resname is None:
                unknown.append(aa)
                continue
            residues.append(resname)
        if unknown:
            report.append(
                f"FASTA record {record_id!r} ignored unsupported residue code(s): {', '.join(sorted(set(unknown)))}."
            )
        if residues:
            sequences.append(Sequence(chain_id, residues))
            report.append(
                f"FASTA sequence loaded for chain {chain_id or '_'} with {len(residues)} residue(s)."
            )
    return sequences, report


def process_structure(
    structure_path: Path,
    options: ProcessingOptions,
    fasta_path: Path | None = None,
) -> ProcessingResult:
    PDBFixer, PDBFile, _PDBxFile, Sequence = _import_openmm_stack()

    report: list[str] = []
    warnings: list[str] = []
    if options.fill_missing_residues and fasta_path is None and structure_path.suffix.lower() == ".pdb":
        try:
            has_seqres = any(
                line.startswith("SEQRES")
                for line in structure_path.read_text(errors="replace").splitlines()
            )
        except OSError:
            has_seqres = True
        if not has_seqres:
            warnings.append(
                "Input PDB has no SEQRES records. PDBFixer may be unable to detect/fill internal missing residues; "
                "use the original PDB with SEQRES records or mmCIF when possible."
            )
    fixer = PDBFixer(filename=str(structure_path))

    if fasta_path is not None:
        fasta_sequences, fasta_report = _load_fasta_sequences(fasta_path, fixer, Sequence)
        report.extend(fasta_report)
        if fasta_sequences:
            fixer.sequences = fasta_sequences
            report.append("PDBFixer sequence metadata was replaced with FASTA sequence input.")

    deferred_mutations, mutation_report = _apply_pdbfixer_mutations(fixer, options.mutations)
    report.extend(mutation_report)

    fixer.findMissingResidues()
    report.extend(_remove_terminal_missing_residue_entries(fixer))
    gaps = _missing_gaps_from_fixer(fixer)
    long_gaps = _long_gaps(gaps, options.long_gap_threshold)
    if gaps:
        for gap in gaps:
            report.append(
                f"Missing residues: chain {gap.chain_id or '_'} after residue {gap.after_residue or '?'} "
                f"has {gap.length} missing residue(s): {', '.join(gap.residue_names) or 'unknown'}."
            )
    else:
        report.append("No missing internal residues reported by PDBFixer.")

    if long_gaps and options.fill_missing_residues and options.long_gap_action == "ask":
        raise LongGapDecisionRequired(long_gaps)

    if not options.fill_missing_residues:
        fixer.missingResidues = {}
        report.append("Missing residue filling disabled by user option.")
    elif long_gaps and options.long_gap_action == "split":
        _remove_long_missing_residue_entries(fixer, options.long_gap_threshold)
        warnings.append(
            "Long missing residue segment(s) were not filled. Output keeps chain IDs and adds SEGID labels around split pieces."
        )
    elif long_gaps and options.long_gap_action == "continue":
        warnings.append("User approved filling missing residue segment(s) longer than threshold.")

    if options.add_missing_atoms:
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        report.append("PDBFixer added missing atoms for known residues.")
    else:
        report.append("Missing atom filling disabled by user option.")

    if options.add_hydrogens:
        fixer.addMissingHydrogens(options.ph)
        report.append(f"PDBFixer added missing hydrogens at pH {options.ph:.2f}.")

    pdb_text = _write_fixer_to_pdb_text(fixer, PDBFile)
    if structure_path.suffix.lower() in {".cif", ".mmcif"}:
        pdb_text, chain_id_report = restore_chain_ids_from_mmcif(pdb_text, structure_path)
        report.extend(chain_id_report)

    if options.use_propka:
        with tempfile.TemporaryDirectory() as tmpdir:
            propka_input = Path(tmpdir) / "propka_input.pdb"
            propka_input.write_text(pdb_text)
            propka_mutations, propka_report = run_propka_prediction(propka_input, options.ph)
            report.extend(propka_report)
            if options.propka_apply_predictions:
                deferred_mutations = tuple(list(deferred_mutations) + propka_mutations)
            elif propka_mutations:
                report.append("PropKa suggestions were reported only; automatic residue-name changes were disabled.")

    pdb_text, zinc_report = apply_zinc_coordination_names(
        pdb_text,
        cutoff_angstrom=options.zn_cutoff_angstrom,
        histidine_mode=options.histidine_mode,
    )
    report.extend(zinc_report)

    if options.amber_first_phosphate_cleanup:
        pdb_text, amber_report, amber_warnings = remove_first_nucleic_phosphate_atoms(pdb_text)
        report.extend(amber_report)
        warnings.extend(amber_warnings)

    if deferred_mutations:
        pdb_text, manual_report = apply_manual_residue_names(pdb_text, deferred_mutations)
        report.extend(manual_report)
        for mutation in deferred_mutations:
            if mutation.new_resname.strip().upper() in STANDARD_AA:
                warnings.append(
                    f"Mutation chain {mutation.chain_id} residue {mutation.residue_number} -> "
                    f"{mutation.new_resname.upper()} was applied as a residue-name replacement only."
                )

    if long_gaps and options.long_gap_action == "split":
        pdb_text, split_report = apply_split_segids(pdb_text, long_gaps)
        report.extend(split_report)
    else:
        report.append("Preserved SEGID fields while preserving chain IDs.")

    if options.add_terminal_caps:
        pdb_text, capping_report, capping_warnings = add_terminal_caps(pdb_text)
        report.extend(capping_report)
        warnings.extend(capping_warnings)
    else:
        report.append("Terminal ACE/NME capping disabled by user option.")

    if options.add_gromacs_ter_records:
        pdb_text, ter_report = insert_gromacs_ter_records(pdb_text)
        report.extend(ter_report)
    else:
        report.append("GROMACS TER separation disabled by user option.")

    output_format = options.output_format.lower()
    if output_format == "cif":
        try:
            output_text = pdb_to_cif_text(pdb_text)
            report.append("Converted final PDB text to mmCIF with Biopython.")
        except Exception as exc:
            output_text = pdb_text
            output_format = "pdb"
            warnings.append(f"CIF export failed, returned PDB instead: {exc}")
    else:
        output_text = pdb_text

    return ProcessingResult(
        output_text=output_text,
        output_format=output_format,
        warnings=warnings,
        report=report,
    )
