from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

from .models import MutationRequest, ProcessingOptions
from .pdb_text import apply_altloc_selections
from .workflow import LongGapDecisionRequired, process_structure


def _load_config(path: Path | None) -> dict:
    if path is None:
        return {}
    return json.loads(path.read_text())


def _mutation_from_dict(raw: dict) -> MutationRequest:
    return MutationRequest(
        chain_id=str(raw.get("chain_id", "")).strip(),
        residue_number=str(raw.get("residue_number", "")).strip(),
        insertion_code=str(raw.get("insertion_code", "")).strip(),
        new_resname=str(raw.get("new_resname", "")).strip().upper(),
    )


def _options_from_config(config: dict) -> ProcessingOptions:
    mutations = tuple(_mutation_from_dict(item) for item in config.get("mutations", []))
    return ProcessingOptions(
        fill_missing_residues=bool(config.get("fill_missing_residues", True)),
        add_missing_atoms=bool(config.get("add_missing_atoms", True)),
        add_hydrogens=bool(config.get("add_hydrogens", False)),
        ph=float(config.get("ph", 7.4)),
        long_gap_threshold=int(config.get("long_gap_threshold", 20)),
        long_gap_action=str(config.get("long_gap_action", "ask")),
        zn_cutoff_angstrom=float(config.get("zn_cutoff_angstrom", 2.8)),
        histidine_mode=str(config.get("histidine_mode", "auto")),
        amber_first_phosphate_cleanup=bool(config.get("amber_first_phosphate_cleanup", True)),
        add_terminal_caps=bool(config.get("add_terminal_caps", True)),
        add_gromacs_ter_records=bool(config.get("add_gromacs_ter_records", True)),
        use_propka=bool(config.get("use_propka", False)),
        propka_apply_predictions=bool(config.get("propka_apply_predictions", False)),
        mutations=mutations,
        output_format=str(config.get("output_format", "pdb")),
    )


def _write_report(path: Path | None, payload: dict) -> None:
    if path is None:
        print(json.dumps(payload, indent=2))
    else:
        path.write_text(json.dumps(payload, indent=2))


def run(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the protein preprocessing workflow.")
    parser.add_argument("--input", required=True, type=Path, help="Input PDB/CIF path")
    parser.add_argument("--output", required=True, type=Path, help="Processed output path")
    parser.add_argument("--config", type=Path, help="JSON config path")
    parser.add_argument("--report", type=Path, help="JSON report output path")
    parser.add_argument("--fasta", type=Path, help="Optional FASTA path")
    args = parser.parse_args(argv)

    try:
        config = _load_config(args.config)
        structure_path = args.input
        pre_report: list[str] = []
        if config.get("apply_altloc_selection", False) and structure_path.suffix.lower() == ".pdb":
            pdb_text = structure_path.read_text(errors="replace")
            selected_text, altloc_report = apply_altloc_selections(
                pdb_text,
                config.get("altloc_selections", {}),
            )
            altloc_path = args.output.with_suffix(".altloc_input.pdb")
            altloc_path.write_text(selected_text)
            structure_path = altloc_path
            pre_report.extend(altloc_report)

        options = _options_from_config(config)
        result = process_structure(structure_path, options, fasta_path=args.fasta)
        args.output.write_text(result.output_text)
        payload = {
            "status": "ok",
            "output": str(args.output),
            "output_format": result.output_format,
            "warnings": result.warnings,
            "report": pre_report + result.report,
        }
        _write_report(args.report, payload)
        return 0
    except LongGapDecisionRequired as exc:
        payload = {
            "status": "long_gap_required",
            "message": str(exc),
            "gaps": [asdict(gap) for gap in exc.gaps],
        }
        _write_report(args.report, payload)
        return 2
    except Exception as exc:
        payload = {
            "status": "error",
            "message": str(exc),
            "exception_type": exc.__class__.__name__,
        }
        _write_report(args.report, payload)
        return 1


def main() -> None:
    sys.exit(run())


if __name__ == "__main__":
    main()
