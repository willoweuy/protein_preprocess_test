from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class MutationRequest:
    chain_id: str
    residue_number: str
    new_resname: str
    insertion_code: str = ""


@dataclass(frozen=True)
class MissingGap:
    chain_id: str
    after_residue: str
    length: int
    residue_names: tuple[str, ...] = ()


@dataclass(frozen=True)
class AltlocResidue:
    key: str
    chain_id: str
    residue_number: str
    insertion_code: str
    resname: str
    choices: tuple[str, ...]
    default_choice: str
    occupancy_by_choice: dict[str, float]


@dataclass(frozen=True)
class ProcessingOptions:
    fill_missing_residues: bool = True
    add_missing_atoms: bool = True
    add_hydrogens: bool = False
    ph: float = 7.4
    long_gap_threshold: int = 20
    long_gap_action: str = "ask"
    zn_cutoff_angstrom: float = 2.8
    histidine_mode: str = "auto"
    amber_first_phosphate_cleanup: bool = True
    add_terminal_caps: bool = True
    add_gromacs_ter_records: bool = True
    use_propka: bool = False
    propka_apply_predictions: bool = False
    mutations: tuple[MutationRequest, ...] = ()
    output_format: str = "pdb"


@dataclass
class ProcessingResult:
    output_text: str
    output_format: str
    warnings: list[str] = field(default_factory=list)
    report: list[str] = field(default_factory=list)
    output_path: Path | None = None
