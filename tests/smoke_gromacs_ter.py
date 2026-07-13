from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_process.pdb_text import insert_gromacs_ter_records


def _atom(serial: int, atom_name: str, resname: str, chain: str, resseq: int, record: str = "ATOM") -> str:
    element = "".join(ch for ch in atom_name if ch.isalpha())[:2].upper()
    return (
        f"{record:<6}{serial:5d} {atom_name:>4} {resname:>3} {chain}{resseq:4d}    "
        f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2}\n"
    )


def main() -> None:
    ligand_ion_chain = "".join(
        [
            _atom(1, "P", "ATP", "A", 1801, "HETATM"),
            _atom(2, "MG", "MG", "A", 1802, "HETATM"),
            _atom(3, "ZN", "ZN", "A", 1805, "HETATM"),
        ]
    )
    separated, report = insert_gromacs_ter_records(ligand_ion_chain)
    assert separated.count("\nTER") == 2, separated
    assert "inserted 2 TER" in report[0]

    protein_chain = "".join(
        [
            _atom(1, "N", "ALA", "A", 1),
            _atom(2, "CA", "ALA", "A", 1),
            _atom(3, "N", "SER", "A", 2),
            _atom(4, "CA", "SER", "A", 2),
        ]
    )
    continuous, report2 = insert_gromacs_ter_records(protein_chain)
    assert "TER" not in continuous
    assert "no additional TER" in report2[0]
    print("smoke gromacs ter ok")


if __name__ == "__main__":
    main()
