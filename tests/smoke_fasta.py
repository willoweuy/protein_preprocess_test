from pathlib import Path
import sys
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_process import ProcessingOptions, process_structure


def main() -> None:
    pdb = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       1.500   1.000   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.500   2.000   0.000  1.00  0.00           O
ATOM      5  N   SER A   3       2.500   1.000   0.000  1.00  0.00           N
ATOM      6  CA  SER A   3       3.500   1.000   0.000  1.00  0.00           C
ATOM      7  C   SER A   3       4.000   2.000   0.000  1.00  0.00           C
ATOM      8  O   SER A   3       4.000   3.000   0.000  1.00  0.00           O
TER
END
"""
    with TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        pdb_path = tmp / "input.pdb"
        fasta_path = tmp / "seq.fasta"
        pdb_path.write_text(pdb)
        fasta_path.write_text(">A\nAGS\n")
        result = process_structure(
            pdb_path,
            ProcessingOptions(
                fill_missing_residues=True,
                add_missing_atoms=False,
                amber_first_phosphate_cleanup=False,
            ),
            fasta_path=fasta_path,
        )
    report_text = "\n".join(result.report)
    assert "FASTA sequence loaded for chain A with 3 residue(s)." in report_text
    assert "has 1 missing residue(s): GLY" in report_text
    print("smoke fasta ok")


if __name__ == "__main__":
    main()
