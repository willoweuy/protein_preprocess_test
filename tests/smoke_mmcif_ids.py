from pathlib import Path
import sys
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_process.pdb_text import get_segid, restore_chain_ids_from_mmcif


def _atom_line_with_chain(chain_id: str) -> str:
    return (
        f"ATOM  {1:5d}  CA  PRO {chain_id}{248:4d}    "
        f"{203.614:8.3f}{175.173:8.3f}{179.353:8.3f}  1.00 179.23           C\n"
    )


def main() -> None:
    mmcif = """data_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1649 C CD . PRO A 1 248 ? 203.614 175.173 179.353 1.00 179.23 ? 248 PRO A CD 1
#
"""
    with TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "input.cif"
        path.write_text(mmcif)
        restored, report = restore_chain_ids_from_mmcif(_atom_line_with_chain("C"), path)

    line = restored.splitlines()[0].ljust(80)
    assert line[21] == "A", restored
    assert get_segid(restored.splitlines(keepends=True)[0]) == "A", restored
    assert "output chain C -> chain A, SEGID A" in report[0], report
    print("smoke mmcif ids ok")


if __name__ == "__main__":
    main()
