from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_process.pdb_text import remove_first_nucleic_phosphate_atoms


def _atom(serial: int, name: str, resname: str = "DA", chain: str = "A", resseq: int = 1) -> str:
    return (
        f"ATOM  {serial:5d} {name:<4} {resname:>3} {chain}{resseq:4d}    "
        f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           {name[0]:>2}\n"
    )


def main() -> None:
    pdb = "".join(
        [
            _atom(1, "P"),
            _atom(2, "O1P"),
            _atom(3, "O2P"),
            _atom(4, "OP1"),
            _atom(5, "OP2"),
            _atom(6, "C1'"),
            _atom(7, "P", resseq=2),
        ]
    )
    cleaned, report, warnings = remove_first_nucleic_phosphate_atoms(pdb)
    assert " C1'" in cleaned
    assert " DA A   2" in cleaned
    for name in (" P  ", "O1P", "O2P", "OP1", "OP2"):
        assert cleaned.count(name) == (1 if name == " P  " else 0), cleaned
    assert "removed O1P, O2P, OP1, OP2, P" in report[0]
    assert warnings and "removed 5 first-phosphate atoms" in warnings[0]
    print("smoke amber cleanup ok")


if __name__ == "__main__":
    main()
