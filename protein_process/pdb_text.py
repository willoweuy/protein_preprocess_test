from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
import math
from pathlib import Path

from .models import AltlocResidue, MissingGap, MutationRequest

ATOM_RECORDS = ("ATOM  ", "HETATM")
NUCLEIC_RESNAMES = {
    "A",
    "C",
    "G",
    "U",
    "DA",
    "DC",
    "DG",
    "DT",
    "DI",
    "RA",
    "RC",
    "RG",
    "RU",
}
WATER_RESNAMES = {"HOH", "WAT", "TIP3", "SOL"}
ION_RESNAMES = {
    "AL",
    "BA",
    "BR",
    "CA",
    "CD",
    "CL",
    "CO",
    "CS",
    "CU",
    "FE",
    "HG",
    "IOD",
    "K",
    "LI",
    "MG",
    "MN",
    "NA",
    "NI",
    "RB",
    "SR",
    "ZN",
}
STANDARD_AA = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
PROTONATION_RESNAMES = {"ASH", "CYM", "GLH", "HID", "HIE", "HIP", "LYN"}
PROTEIN_RESNAMES = STANDARD_AA | PROTONATION_RESNAMES


@dataclass(frozen=True)
class AtomLine:
    index: int
    line: str
    record: str
    atom_name: str
    resname: str
    chain_id: str
    resseq: str
    insertion_code: str
    x: float
    y: float
    z: float
    element: str

    @property
    def residue_key(self) -> tuple[str, str, str]:
        return (self.chain_id, self.resseq, self.insertion_code)


@dataclass(frozen=True)
class CifChainIds:
    label_asym_id: str
    auth_asym_id: str


def _pad_pdb_line(line: str) -> str:
    return line.rstrip("\n").ljust(80)


def is_atom_line(line: str) -> bool:
    return line.startswith(ATOM_RECORDS)


def parse_atom_line(index: int, line: str) -> AtomLine | None:
    if not is_atom_line(line):
        return None
    padded = _pad_pdb_line(line)
    try:
        x = float(padded[30:38])
        y = float(padded[38:46])
        z = float(padded[46:54])
    except ValueError:
        return None
    atom_name = padded[12:16].strip()
    element = padded[76:78].strip() or "".join(ch for ch in atom_name if ch.isalpha())[:2].upper()
    return AtomLine(
        index=index,
        line=line,
        record=padded[0:6],
        atom_name=atom_name,
        resname=padded[17:20].strip(),
        chain_id=padded[21].strip(),
        resseq=padded[22:26].strip(),
        insertion_code=padded[26].strip(),
        x=x,
        y=y,
        z=z,
        element=element.upper(),
    )


def iter_atom_lines(lines: list[str]) -> list[AtomLine]:
    atoms: list[AtomLine] = []
    for index, line in enumerate(lines):
        atom = parse_atom_line(index, line)
        if atom is not None:
            atoms.append(atom)
    return atoms


def _altloc_from_line(line: str) -> str:
    if not is_atom_line(line):
        return ""
    return _pad_pdb_line(line)[16].strip()


def _occupancy_from_line(line: str) -> float:
    if not is_atom_line(line):
        return 0.0
    try:
        return float(_pad_pdb_line(line)[54:60])
    except ValueError:
        return 0.0


def _residue_key_string(chain_id: str, resseq: str, insertion_code: str, resname: str) -> str:
    return f"{chain_id}|{resseq}|{insertion_code}|{resname}"


def _split_residue_key_string(key: str) -> tuple[str, str, str, str]:
    chain_id, resseq, insertion_code, resname = (key.split("|") + ["", "", "", ""])[:4]
    return chain_id, resseq, insertion_code, resname


def find_altloc_residues(pdb_text: str) -> tuple[AltlocResidue, ...]:
    lines = pdb_text.splitlines(keepends=True)
    occupancy_values: dict[str, dict[str, list[float]]] = {}
    residue_labels: dict[str, tuple[str, str, str, str]] = {}

    for index, line in enumerate(lines):
        atom = parse_atom_line(index, line)
        if atom is None:
            continue
        altloc = _altloc_from_line(line)
        if not altloc:
            continue
        key = _residue_key_string(atom.chain_id, atom.resseq, atom.insertion_code, atom.resname)
        occupancy_values.setdefault(key, {}).setdefault(altloc, []).append(_occupancy_from_line(line))
        residue_labels[key] = (atom.chain_id, atom.resseq, atom.insertion_code, atom.resname)

    residues: list[AltlocResidue] = []
    for key, by_altloc in sorted(occupancy_values.items()):
        if len(by_altloc) < 2:
            continue
        occupancy_by_choice = {
            altloc: sum(values) / len(values) if values else 0.0
            for altloc, values in by_altloc.items()
        }
        max_occupancy = max(occupancy_by_choice.values())
        tied = [
            altloc
            for altloc, occupancy in occupancy_by_choice.items()
            if abs(occupancy - max_occupancy) < 1e-6
        ]
        default_choice = "A" if "A" in tied else sorted(tied)[0]
        chain_id, resseq, insertion_code, resname = residue_labels[key]
        residues.append(
            AltlocResidue(
                key=key,
                chain_id=chain_id,
                residue_number=resseq,
                insertion_code=insertion_code,
                resname=resname,
                choices=tuple(sorted(occupancy_by_choice)),
                default_choice=default_choice,
                occupancy_by_choice=occupancy_by_choice,
            )
        )
    return tuple(residues)


def _clear_altloc(line: str) -> str:
    padded = _pad_pdb_line(line)
    fixed = f"{padded[:16]} {padded[17:]}"
    return fixed.rstrip() + "\n"


def apply_altloc_selections(
    pdb_text: str,
    selections: dict[str, str] | None = None,
) -> tuple[str, list[str]]:
    altloc_residues = find_altloc_residues(pdb_text)
    if not altloc_residues:
        return pdb_text, ["Alternate occupancy: no PDB alternate-location records detected."]

    selected_by_key = {
        residue.key: (selections or {}).get(residue.key, residue.default_choice)
        for residue in altloc_residues
    }
    altloc_keys = set(selected_by_key)
    lines = pdb_text.splitlines(keepends=True)
    kept_lines: list[str] = []

    for index, line in enumerate(lines):
        atom = parse_atom_line(index, line)
        if atom is None:
            kept_lines.append(line)
            continue
        key = _residue_key_string(atom.chain_id, atom.resseq, atom.insertion_code, atom.resname)
        altloc = _altloc_from_line(line)
        if key not in altloc_keys or not altloc:
            kept_lines.append(line)
            continue
        if altloc == selected_by_key[key]:
            kept_lines.append(_clear_altloc(line))

    report: list[str] = []
    by_key = {residue.key: residue for residue in altloc_residues}
    for key, selected in selected_by_key.items():
        residue = by_key[key]
        occupancy = residue.occupancy_by_choice.get(selected, 0.0)
        options = ", ".join(
            f"{choice}:{residue.occupancy_by_choice.get(choice, 0.0):.2f}"
            for choice in residue.choices
        )
        report.append(
            f"Alternate occupancy: chain {residue.chain_id or '_'} residue "
            f"{residue.residue_number}{residue.insertion_code} {residue.resname} selected "
            f"{selected} (mean occupancy {occupancy:.2f}; options {options})."
        )
    return "".join(kept_lines), report


def set_resname(line: str, resname: str) -> str:
    padded = _pad_pdb_line(line)
    fixed = f"{padded[:17]}{resname.upper():>3}{padded[20:]}"
    return fixed.rstrip() + "\n"


def set_segid(line: str, segid: str) -> str:
    padded = _pad_pdb_line(line)
    fixed = f"{padded[:72]}{segid[:4]:<4}{padded[76:]}"
    return fixed.rstrip() + "\n"


def set_chain_id(line: str, chain_id: str) -> str:
    if not is_atom_line(line):
        return line
    padded = _pad_pdb_line(line)
    pdb_chain_id = (chain_id.strip() or " ")[:1]
    fixed = f"{padded[:21]}{pdb_chain_id}{padded[22:]}"
    return fixed.rstrip() + "\n"


def _set_serial(line: str, serial: int) -> str:
    padded = _pad_pdb_line(line)
    fixed = f"{padded[:6]}{serial:5d}{padded[11:]}"
    return fixed.rstrip() + "\n"


def get_segid(line: str) -> str:
    if not is_atom_line(line):
        return ""
    return _pad_pdb_line(line)[72:76].strip()


def _serial_from_line(line: str) -> int:
    if not is_atom_line(line):
        return 0
    try:
        return int(_pad_pdb_line(line)[6:11])
    except ValueError:
        return 0


def _format_pdb_atom_line(
    serial: int,
    atom_name: str,
    resname: str,
    chain_id: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    segid: str,
    element: str,
) -> str:
    line = (
        f"ATOM  {serial:5d} {atom_name:>4} {resname:>3} {(chain_id or ' ')[:1]}"
        f"{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00           {element:>2}\n"
    )
    return set_segid(line, segid)


def _format_ter_line(serial: int, atom: AtomLine) -> str:
    return (
        f"TER   {serial:5d}      {atom.resname:>3} {(atom.chain_id or ' ')[:1]}"
        f"{_residue_number_or_default(atom.resseq, 1):4d}{(atom.insertion_code or ' ')[:1]}\n"
    )


def _v_add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _v_sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _v_scale(a: tuple[float, float, float], scale: float) -> tuple[float, float, float]:
    return (a[0] * scale, a[1] * scale, a[2] * scale)


def _v_norm(a: tuple[float, float, float]) -> tuple[float, float, float]:
    length = math.sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)
    if length < 1e-8:
        return (1.0, 0.0, 0.0)
    return (a[0] / length, a[1] / length, a[2] / length)


def _v_cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _v_dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _v_length(a: tuple[float, float, float]) -> float:
    return math.sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)


def _perpendicular(a: tuple[float, float, float]) -> tuple[float, float, float]:
    perp = _v_cross(a, (0.0, 0.0, 1.0))
    if _v_length(perp) < 1e-8:
        perp = _v_cross(a, (0.0, 1.0, 0.0))
    return _v_norm(perp)


def _plane_perpendicular_axis(
    bond_axis: tuple[float, float, float],
    point_a: tuple[float, float, float],
    point_b: tuple[float, float, float],
    point_c: tuple[float, float, float] | None,
) -> tuple[float, float, float]:
    if point_c is not None:
        normal = _v_cross(_v_sub(point_b, point_a), _v_sub(point_c, point_a))
        if _v_length(normal) >= 1e-8:
            in_plane = _v_cross(_v_norm(normal), bond_axis)
            if _v_length(in_plane) >= 1e-8:
                return _v_norm(in_plane)
    return _perpendicular(bond_axis)


def _trigonal_directions(
    bond_axis: tuple[float, float, float],
    in_plane_axis: tuple[float, float, float],
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    bond_axis = _v_norm(bond_axis)
    in_plane_axis = _perpendicular_component(bond_axis, in_plane_axis)
    root_three_over_two = math.sqrt(3.0) / 2.0
    opposite = _v_scale(bond_axis, -0.5)
    side = _v_scale(in_plane_axis, root_three_over_two)
    return _v_norm(_v_add(opposite, side)), _v_norm(_v_sub(opposite, side))


def _perpendicular_component(
    axis: tuple[float, float, float],
    reference: tuple[float, float, float],
) -> tuple[float, float, float]:
    axis = _v_norm(axis)
    component = _v_sub(reference, _v_scale(axis, _v_dot(reference, axis)))
    if _v_length(component) < 1e-8:
        return _perpendicular(axis)
    return _v_norm(component)


def _direction_120_from(
    axis: tuple[float, float, float],
    in_plane_axis: tuple[float, float, float],
) -> tuple[float, float, float]:
    axis = _v_norm(axis)
    in_plane_axis = _perpendicular_component(axis, in_plane_axis)
    return _v_norm(
        _v_add(
            _v_scale(axis, -0.5),
            _v_scale(in_plane_axis, math.sqrt(3.0) / 2.0),
        )
    )


def _atom_position(atom: AtomLine) -> tuple[float, float, float]:
    return (atom.x, atom.y, atom.z)


def _residue_sort_key(key: tuple[str, str, str]) -> tuple[int, str, str]:
    _chain_id, resseq, insertion_code = key
    try:
        return (int(resseq), insertion_code, resseq)
    except ValueError:
        return (0, insertion_code, resseq)


def _residue_number_or_default(resseq: str, default: int) -> int:
    try:
        return int(resseq)
    except ValueError:
        return default


def _gromacs_residue_type(resname: str) -> str:
    name = resname.upper()
    if name in PROTEIN_RESNAMES or name in {"ACE", "NMA", "NME"}:
        return "protein"
    if name in NUCLEIC_RESNAMES:
        return "nucleic"
    if name in ION_RESNAMES:
        return "ion"
    if name in WATER_RESNAMES:
        return "water"
    return "ligand"


def _needs_gromacs_ter(previous: AtomLine, current: AtomLine) -> bool:
    if previous.chain_id != current.chain_id:
        return False
    prev_type = _gromacs_residue_type(previous.resname)
    curr_type = _gromacs_residue_type(current.resname)
    non_polymer_types = {"ligand", "ion", "water"}
    if prev_type in non_polymer_types or curr_type in non_polymer_types:
        return True
    return prev_type != curr_type


def clear_segid(line: str) -> str:
    if not is_atom_line(line):
        return line
    padded = _pad_pdb_line(line)
    fixed = f"{padded[:72]}    {padded[76:]}"
    return fixed.rstrip() + "\n"


def remove_all_segids(pdb_text: str) -> str:
    return "".join(clear_segid(line) for line in pdb_text.splitlines(keepends=True))


def _as_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def _clean_cif_id(value: str) -> str:
    value = value.strip()
    if value in {"", ".", "?"}:
        return ""
    return value


def read_mmcif_chain_ids(mmcif_path: Path) -> tuple[CifChainIds, ...]:
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    except ImportError as exc:
        raise RuntimeError("Biopython is required for mmCIF chain ID restoration.") from exc

    mmcif = MMCIF2Dict(str(mmcif_path))
    label_ids = _as_list(mmcif.get("_atom_site.label_asym_id"))
    auth_ids = _as_list(mmcif.get("_atom_site.auth_asym_id"))
    if not label_ids or not auth_ids:
        return ()

    chain_ids: list[CifChainIds] = []
    seen: set[tuple[str, str]] = set()
    for label_id, auth_id in zip(label_ids, auth_ids):
        label = _clean_cif_id(label_id)
        auth = _clean_cif_id(auth_id)
        if not label and not auth:
            continue
        key = (label, auth)
        if key in seen:
            continue
        seen.add(key)
        chain_ids.append(CifChainIds(label_asym_id=label, auth_asym_id=auth))
    return tuple(chain_ids)


def _pdb_chain_order(lines: list[str]) -> list[str]:
    order: list[str] = []
    seen: set[str] = set()
    for atom in iter_atom_lines(lines):
        if atom.chain_id in seen:
            continue
        seen.add(atom.chain_id)
        order.append(atom.chain_id)
    return order


def restore_chain_ids_from_mmcif(pdb_text: str, mmcif_path: Path) -> tuple[str, list[str]]:
    chain_ids = read_mmcif_chain_ids(mmcif_path)
    if not chain_ids:
        return pdb_text, ["mmCIF chain ID restoration: no label/auth asym IDs found."]

    lines = pdb_text.splitlines(keepends=True)
    pdb_chain_order = _pdb_chain_order(lines)
    if not pdb_chain_order:
        return pdb_text, ["mmCIF chain ID restoration: no PDB ATOM/HETATM chains found."]

    records_by_id: dict[str, CifChainIds] = {}
    for record in chain_ids:
        if record.label_asym_id:
            records_by_id.setdefault(record.label_asym_id, record)
        if record.auth_asym_id:
            records_by_id.setdefault(record.auth_asym_id, record)

    used_records: set[CifChainIds] = set()
    mapping: dict[str, CifChainIds] = {}
    remaining_records = list(chain_ids)
    for current_chain_id in pdb_chain_order:
        record = records_by_id.get(current_chain_id)
        if record in used_records:
            record = None
        if record is None:
            while remaining_records and remaining_records[0] in used_records:
                remaining_records.pop(0)
            record = remaining_records[0] if remaining_records else None
        if record is None:
            continue
        mapping[current_chain_id] = record
        used_records.add(record)

    if not mapping:
        return pdb_text, ["mmCIF chain ID restoration: no chain mapping could be inferred."]

    report: list[str] = []
    truncations: set[str] = set()
    for index, line in enumerate(lines):
        atom = parse_atom_line(index, line)
        if atom is None:
            continue
        record = mapping.get(atom.chain_id)
        if record is None:
            continue
        chain_id = record.auth_asym_id or record.label_asym_id or atom.chain_id
        segid = record.label_asym_id or record.auth_asym_id or get_segid(line)
        if len(chain_id.strip()) > 1:
            truncations.add(chain_id)
        updated = set_chain_id(line, chain_id)
        updated = set_segid(updated, segid)
        lines[index] = updated

    for current_chain_id, record in mapping.items():
        chain_id = record.auth_asym_id or record.label_asym_id or current_chain_id
        segid = record.label_asym_id or record.auth_asym_id or ""
        report.append(
            f"mmCIF chain ID restoration: output chain {current_chain_id or '_'} -> "
            f"chain {chain_id[:1] or '_'}, SEGID {segid[:4] or '_'} "
            f"(auth_asym_id={record.auth_asym_id or '_'}, label_asym_id={record.label_asym_id or '_'})."
        )
    for chain_id in sorted(truncations):
        report.append(
            f"mmCIF chain ID restoration: auth_asym_id {chain_id!r} was truncated to {chain_id[:1]!r} for PDB output."
        )
    return "".join(lines), report


def _split_number_segid(chain_id: str, segment_index: int) -> str:
    chain_prefix = (chain_id.strip() or "X")[:1]
    if segment_index < 10:
        return f"{chain_prefix}{segment_index}"
    return f"{chain_prefix}{segment_index:03d}"[-4:]


def squared_distance(a: AtomLine, b: AtomLine) -> float:
    return (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2


def apply_zinc_coordination_names(
    pdb_text: str,
    cutoff_angstrom: float,
    histidine_mode: str,
) -> tuple[str, list[str]]:
    lines = pdb_text.splitlines(keepends=True)
    atoms = iter_atom_lines(lines)
    zn_atoms = [
        atom
        for atom in atoms
        if atom.resname.upper() in {"ZN", "ZN2"} or atom.element.upper() == "ZN"
    ]
    if not zn_atoms:
        return pdb_text, ["No Zn atoms detected."]

    cutoff2 = cutoff_angstrom * cutoff_angstrom
    residue_updates: dict[tuple[str, str, str], str] = {}
    his_sites: dict[tuple[str, str, str], set[str]] = {}
    report: list[str] = [f"Detected {len(zn_atoms)} Zn atom(s)."]

    for zn in zn_atoms:
        for atom in atoms:
            if atom.index == zn.index:
                continue
            if squared_distance(zn, atom) > cutoff2:
                continue
            resname = atom.resname.upper()
            atom_name = atom.atom_name.upper()
            if resname == "CYS" and (atom_name == "SG" or atom.element == "S"):
                residue_updates[atom.residue_key] = "CYM"
                report.append(
                    f"Zn coordination: chain {atom.chain_id or '_'} residue {atom.resseq}{atom.insertion_code} CYS -> CYM."
                )
            elif resname in {"HIS", "HID", "HIE", "HIP"} and atom_name in {"ND1", "NE2"}:
                his_sites.setdefault(atom.residue_key, set()).add(atom_name)

    for key, donor_atoms in his_sites.items():
        if histidine_mode in {"HID", "HIE", "HIP"}:
            new_name = histidine_mode
        elif donor_atoms == {"ND1"}:
            new_name = "HID"
        elif donor_atoms == {"NE2"}:
            new_name = "HIE"
        else:
            new_name = "HIP"
        residue_updates[key] = new_name
        chain_id, resseq, insertion_code = key
        donors = ",".join(sorted(donor_atoms))
        report.append(
            f"Zn coordination: chain {chain_id or '_'} residue {resseq}{insertion_code} HIS donors {donors} -> {new_name}."
        )

    if not residue_updates:
        report.append(f"No CYS/HIS Zn coordination found within {cutoff_angstrom:.2f} A.")
        return pdb_text, report

    for atom in atoms:
        new_resname = residue_updates.get(atom.residue_key)
        if new_resname:
            lines[atom.index] = set_resname(lines[atom.index], new_resname)
    return "".join(lines), report


def apply_manual_residue_names(
    pdb_text: str,
    mutations: tuple[MutationRequest, ...],
) -> tuple[str, list[str]]:
    if not mutations:
        return pdb_text, []

    lines = pdb_text.splitlines(keepends=True)
    atoms = iter_atom_lines(lines)
    requests = {
        (
            mutation.chain_id.strip(),
            str(mutation.residue_number).strip(),
            mutation.insertion_code.strip(),
        ): mutation.new_resname.strip().upper()
        for mutation in mutations
        if mutation.chain_id.strip() and str(mutation.residue_number).strip() and mutation.new_resname.strip()
    }
    touched: set[tuple[str, str, str]] = set()
    for atom in atoms:
        new_resname = requests.get(atom.residue_key)
        if new_resname:
            lines[atom.index] = set_resname(lines[atom.index], new_resname)
            touched.add(atom.residue_key)

    report: list[str] = []
    for key, new_resname in requests.items():
        chain_id, resseq, insertion_code = key
        if key in touched:
            report.append(f"Manual residue name: chain {chain_id} residue {resseq}{insertion_code} -> {new_resname}.")
        else:
            report.append(f"Manual residue name skipped: chain {chain_id} residue {resseq}{insertion_code} was not found.")
    return "".join(lines), report


def add_terminal_caps(pdb_text: str) -> tuple[str, list[str], list[str]]:
    lines = pdb_text.splitlines(keepends=True)
    atoms = iter_atom_lines(lines)
    if not atoms:
        return pdb_text, ["Terminal capping: no atoms found."], []

    max_serial = max((_serial_from_line(line) for line in lines), default=0)
    residues_by_segment: dict[tuple[str, str], dict[tuple[str, str, str], list[AtomLine]]] = {}
    segment_has_ace: set[tuple[str, str]] = set()
    segment_has_nma: set[tuple[str, str]] = set()
    atom_segids = {atom.index: get_segid(lines[atom.index]) for atom in atoms}

    for atom in atoms:
        segid = atom_segids.get(atom.index, "")
        segment_key = (atom.chain_id, segid)
        resname = atom.resname.upper()
        if resname == "ACE":
            segment_has_ace.add(segment_key)
            continue
        if resname in {"NMA", "NME"}:
            segment_has_nma.add(segment_key)
            continue
        if resname not in PROTEIN_RESNAMES:
            continue
        residues_by_segment.setdefault(segment_key, {}).setdefault(atom.residue_key, []).append(atom)

    if not residues_by_segment:
        return pdb_text, ["Terminal capping: no protein segments found."], []

    insert_before: dict[int, list[str]] = {}
    insert_after: dict[int, list[str]] = {}
    remove_line_indices: set[int] = set()
    report: list[str] = []
    warnings: list[str] = []

    for (chain_id, segid), residues in residues_by_segment.items():
        ordered_keys = sorted(residues, key=_residue_sort_key)
        if not ordered_keys:
            continue

        first_key = ordered_keys[0]
        last_key = ordered_keys[-1]
        first_atoms = residues[first_key]
        last_atoms = residues[last_key]
        first_by_name = {atom.atom_name.upper(): atom for atom in first_atoms}
        last_by_name = {atom.atom_name.upper(): atom for atom in last_atoms}
        segment_label = f"chain {chain_id or '_'} SEGID {segid or '_'}"

        if (chain_id, segid) not in segment_has_ace:
            n_atom = first_by_name.get("N")
            ca_atom = first_by_name.get("CA")
            c_atom = first_by_name.get("C")
            if n_atom is None or ca_atom is None:
                warnings.append(f"Terminal capping: skipped ACE for {segment_label}; missing N or CA atom.")
            else:
                max_serial += 1
                n_pos = _atom_position(n_atom)
                ca_pos = _atom_position(ca_atom)
                n_to_ca = _v_norm(_v_sub(ca_pos, n_pos))
                n_plane_reference = _v_sub(_atom_position(c_atom), n_pos) if c_atom is not None else _perpendicular(n_to_ca)
                n_to_carbon = _direction_120_from(n_to_ca, n_plane_reference)
                carbon = _v_add(n_pos, _v_scale(n_to_carbon, 1.33))
                bond_to_n = _v_norm(_v_sub(n_pos, carbon))
                in_plane = _perpendicular_component(
                    bond_to_n,
                    _v_sub(ca_pos, carbon),
                )
                option_a, option_b = _trigonal_directions(bond_to_n, in_plane)
                if _v_length(_v_sub(_v_add(carbon, _v_scale(option_a, 1.23)), ca_pos)) >= _v_length(
                    _v_sub(_v_add(carbon, _v_scale(option_b, 1.23)), ca_pos)
                ):
                    oxygen_direction, methyl_direction = option_a, option_b
                else:
                    oxygen_direction, methyl_direction = option_b, option_a
                oxygen = _v_add(carbon, _v_scale(oxygen_direction, 1.23))
                methyl = _v_add(carbon, _v_scale(methyl_direction, 1.50))
                ace_resseq = _residue_number_or_default(first_key[1], 1) - 1
                ace_lines = [
                    _format_pdb_atom_line(max_serial, "C", "ACE", chain_id, ace_resseq, *carbon, segid=segid, element="C"),
                    _format_pdb_atom_line(max_serial + 1, "CH3", "ACE", chain_id, ace_resseq, *methyl, segid=segid, element="C"),
                    _format_pdb_atom_line(max_serial + 2, "O", "ACE", chain_id, ace_resseq, *oxygen, segid=segid, element="O"),
                ]
                max_serial += 2
                first_line_index = min(atom.index for atom in first_atoms)
                insert_before.setdefault(first_line_index, []).extend(ace_lines)
                report.append(f"Terminal capping: added ACE before {segment_label} residue {first_key[1]}{first_key[2]}.")
        else:
            report.append(f"Terminal capping: ACE already present for {segment_label}; skipped.")

        if (chain_id, segid) not in segment_has_nma:
            c_atom = last_by_name.get("C")
            ca_atom = last_by_name.get("CA")
            if c_atom is None or ca_atom is None:
                warnings.append(f"Terminal capping: skipped NME for {segment_label}; missing C or CA atom.")
            else:
                oxt_atom = last_by_name.get("OXT")
                if oxt_atom is not None:
                    remove_line_indices.add(oxt_atom.index)
                max_serial += 1
                c_pos = _atom_position(c_atom)
                ca_pos = _atom_position(ca_atom)
                o_atom = last_by_name.get("O")
                bond_to_ca = _v_norm(_v_sub(ca_pos, c_pos))
                in_plane = _plane_perpendicular_axis(
                    bond_to_ca,
                    c_pos,
                    ca_pos,
                    _atom_position(o_atom) if o_atom is not None else None,
                )
                option_a, option_b = _trigonal_directions(bond_to_ca, in_plane)
                if o_atom is not None:
                    o_direction = _v_norm(_v_sub(_atom_position(o_atom), c_pos))
                    nitrogen_direction = option_a if _v_dot(option_a, o_direction) < _v_dot(option_b, o_direction) else option_b
                else:
                    nitrogen_direction = option_a
                nitrogen = _v_add(c_pos, _v_scale(nitrogen_direction, 1.33))
                n_to_c = _v_norm(_v_sub(c_pos, nitrogen))
                methyl_plane = _perpendicular_component(n_to_c, _v_sub(ca_pos, nitrogen))
                methyl_option_a, methyl_option_b = _trigonal_directions(n_to_c, methyl_plane)
                methyl_direction = methyl_option_a if _v_dot(methyl_option_a, bond_to_ca) < _v_dot(methyl_option_b, bond_to_ca) else methyl_option_b
                methyl = _v_add(nitrogen, _v_scale(methyl_direction, 1.45))
                nma_resseq = _residue_number_or_default(last_key[1], 1) + 1
                nma_lines = [
                    _format_pdb_atom_line(max_serial, "N", "NME", chain_id, nma_resseq, *nitrogen, segid=segid, element="N"),
                    _format_pdb_atom_line(max_serial + 1, "CH3", "NME", chain_id, nma_resseq, *methyl, segid=segid, element="C"),
                ]
                max_serial += 1
                last_line_index = max(atom.index for atom in last_atoms)
                insert_after.setdefault(last_line_index, []).extend(nma_lines)
                extra = " and removed OXT" if oxt_atom is not None else ""
                report.append(f"Terminal capping: added NME after {segment_label} residue {last_key[1]}{last_key[2]}{extra}.")
        else:
            report.append(f"Terminal capping: NMA/NME already present for {segment_label}; skipped.")

    capped_lines: list[str] = []
    for index, line in enumerate(lines):
        capped_lines.extend(insert_before.get(index, []))
        if index not in remove_line_indices:
            capped_lines.append(line)
        capped_lines.extend(insert_after.get(index, []))
    return "".join(capped_lines), report, warnings


def insert_gromacs_ter_records(pdb_text: str) -> tuple[str, list[str]]:
    lines = pdb_text.splitlines(keepends=True)
    max_serial = max((_serial_from_line(line) for line in lines), default=0)
    output: list[str] = []
    previous_atom: AtomLine | None = None
    previous_residue_key: tuple[str, str, str] | None = None
    last_output_was_ter = False
    inserted = 0

    for index, line in enumerate(lines):
        atom = parse_atom_line(index, line)
        if line.startswith("TER"):
            output.append(line)
            previous_atom = None
            previous_residue_key = None
            last_output_was_ter = True
            continue
        if atom is None:
            output.append(line)
            last_output_was_ter = False
            continue

        current_residue_key = atom.residue_key
        if (
            previous_atom is not None
            and previous_residue_key is not None
            and current_residue_key != previous_residue_key
            and _needs_gromacs_ter(previous_atom, atom)
            and not last_output_was_ter
        ):
            max_serial += 1
            output.append(_format_ter_line(max_serial, previous_atom))
            inserted += 1

        output.append(line)
        previous_atom = atom
        previous_residue_key = current_residue_key
        last_output_was_ter = False

    if inserted:
        return "".join(output), [f"GROMACS TER separation: inserted {inserted} TER record(s) at molecule-type boundaries."]
    return pdb_text, ["GROMACS TER separation: no additional TER records were needed."]


def remove_first_nucleic_phosphate_atoms(pdb_text: str) -> tuple[str, list[str], list[str]]:
    lines = pdb_text.splitlines(keepends=True)
    atoms = iter_atom_lines(lines)
    nucleic_residues: dict[str, list[tuple[int, str, str]]] = {}
    for atom in atoms:
        if atom.resname.upper() in NUCLEIC_RESNAMES:
            try:
                residue_number = int(atom.resseq)
            except ValueError:
                continue
            nucleic_residues.setdefault(atom.chain_id, []).append(
                (residue_number, atom.resseq, atom.insertion_code)
            )

    first_by_chain: set[tuple[str, str, str]] = set()
    for chain_id, residue_ids in nucleic_residues.items():
        if not residue_ids:
            continue
        residue_number, resseq, insertion_code = min(residue_ids, key=lambda item: (item[0], item[2]))
        first_by_chain.add((chain_id, resseq, insertion_code))

    if not first_by_chain:
        return pdb_text, ["Amber first-phosphate cleanup: no DNA/RNA residues detected."], []

    remove_names = {"P", "O1P", "O2P", "OP1", "OP2"}
    removed: dict[tuple[str, str, str], list[str]] = {}
    keep_lines: list[str] = []
    for index, line in enumerate(lines):
        atom = parse_atom_line(index, line)
        if atom and atom.residue_key in first_by_chain and atom.atom_name.upper() in remove_names:
            removed.setdefault(atom.residue_key, []).append(atom.atom_name.upper())
            continue
        keep_lines.append(line)

    report: list[str] = []
    warnings: list[str] = []
    for chain_id, resseq, insertion_code in sorted(first_by_chain):
        names = sorted(set(removed.get((chain_id, resseq, insertion_code), [])))
        if names:
            report.append(
                f"Amber cleanup: removed {', '.join(names)} from first nucleic residue chain {chain_id or '_'} {resseq}{insertion_code}."
            )
            if len(names) > 3:
                warnings.append(
                    f"Amber cleanup removed {len(names)} first-phosphate atoms "
                    f"({', '.join(names)}) from first nucleic residue chain {chain_id or '_'} {resseq}{insertion_code}."
                )
        else:
            report.append(
                f"Amber cleanup: first nucleic residue chain {chain_id or '_'} {resseq}{insertion_code} had no P/O1P/O2P/OP1/OP2 atoms to remove."
            )
    return "".join(keep_lines), report, warnings


def apply_split_segids(pdb_text: str, long_gaps: tuple[MissingGap, ...]) -> tuple[str, list[str]]:
    if not long_gaps:
        return pdb_text, []

    lines = pdb_text.splitlines(keepends=True)
    atoms = iter_atom_lines(lines)
    cutoffs_by_chain: dict[str, list[int]] = {}
    for gap in long_gaps:
        try:
            cutoff = int(gap.after_residue)
        except ValueError:
            continue
        cutoffs_by_chain.setdefault(gap.chain_id, []).append(cutoff)
    for cutoffs in cutoffs_by_chain.values():
        cutoffs.sort()

    report: list[str] = []
    for atom in atoms:
        cutoffs = cutoffs_by_chain.get(atom.chain_id)
        if not cutoffs:
            continue
        try:
            residue_number = int(atom.resseq)
        except ValueError:
            continue
        segment_index = 1
        for cutoff in cutoffs:
            if residue_number > cutoff:
                segment_index += 1
        segid = _split_number_segid(atom.chain_id, segment_index)
        lines[atom.index] = set_segid(lines[atom.index], segid)

    for chain_id, cutoffs in cutoffs_by_chain.items():
        assigned = [
            _split_number_segid(chain_id, segment_index)
            for segment_index in range(1, len(cutoffs) + 2)
        ]
        report.append(
            f"Split mode: assigned SEGID labels around chain {chain_id or '_'} gap cutoff(s) "
            f"{', '.join(map(str, cutoffs))}: {', '.join(assigned)}."
        )
    return "".join(lines), report


def pdb_to_cif_text(pdb_text: str, structure_id: str = "processed") -> str:
    try:
        from Bio.PDB import MMCIFIO, PDBParser
    except ImportError as exc:
        raise RuntimeError("Biopython is required for CIF export.") from exc

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(structure_id, StringIO(pdb_text))
    output = StringIO()
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output)
    return output.getvalue()
