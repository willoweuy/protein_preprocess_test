from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_process.pdb_text import apply_altloc_selections, find_altloc_residues


def main() -> None:
    pdb = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA AALA A   1       1.000   0.000   0.000  0.50  0.00           C
ATOM      3  CA BALA A   1       2.000   0.000   0.000  0.50  0.00           C
ATOM      4  CB AALA A   1       1.000   1.000   0.000  0.50  0.00           C
ATOM      5  CB BALA A   1       2.000   1.000   0.000  0.50  0.00           C
ATOM      6  C   ALA A   1       1.000   0.000   1.000  1.00  0.00           C
END
"""
    residues = find_altloc_residues(pdb)
    assert len(residues) == 1
    assert residues[0].default_choice == "A"
    selected, report = apply_altloc_selections(pdb, {})
    assert "CA AALA" not in selected
    assert "CA BALA" not in selected
    assert "CA  ALA" in selected
    assert "2.000   0.000   0.000" not in selected
    assert "selected A" in report[0]
    selected_b, _report_b = apply_altloc_selections(pdb, {residues[0].key: "B"})
    assert "1.000   0.000   0.000" not in selected_b
    assert "2.000   0.000   0.000" in selected_b
    print("smoke altloc ok")


if __name__ == "__main__":
    main()
