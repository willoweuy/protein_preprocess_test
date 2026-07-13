from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_process.models import MissingGap
from protein_process.pdb_text import apply_split_segids, get_segid, set_segid


def _atom_line(serial: int, chain_id: str, residue_number: int) -> str:
    return (
        f"ATOM  {serial:5d}  CA  ALA {chain_id}{residue_number:4d}    "
        f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n"
    )


def _atom_segids(pdb_text: str) -> list[str]:
    return [
        get_segid(line)
        for line in pdb_text.splitlines(keepends=True)
        if line.startswith("ATOM")
    ]


def main() -> None:
    gap = (MissingGap(chain_id="D", after_residue="10", length=20),)

    no_old = _atom_line(1, "D", 5) + _atom_line(2, "D", 15)
    split, report = apply_split_segids(no_old, gap)
    assert _atom_segids(split) == ["D1", "D2"], split
    assert "D1, D2" in report[0]

    with_old = set_segid(_atom_line(1, "D", 5), "segF") + set_segid(_atom_line(2, "D", 15), "segF")
    split_old, report_old = apply_split_segids(with_old, gap)
    assert _atom_segids(split_old) == ["D1", "D2"], split_old
    assert "D1, D2" in report_old[0]
    print("smoke segid ok")


if __name__ == "__main__":
    main()
