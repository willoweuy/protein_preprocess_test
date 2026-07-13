from pathlib import Path
import sys
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_process import ProcessingOptions, process_structure


def main() -> None:
    pdb = """ATOM      1  N   HIS A  12      -1.200   0.000   0.000  1.00  0.00           N
ATOM      2  CA  HIS A  12      -0.300   0.800   0.000  1.00  0.00           C
ATOM      3  C   HIS A  12       0.900   0.000   0.000  1.00  0.00           C
ATOM      4  O   HIS A  12       1.900   0.500   0.000  1.00  0.00           O
ATOM      5  NE2 HIS A  12       0.000   2.200   0.000  1.00  0.00           N
ATOM      6  N   CYS A  13       0.800  -1.300   0.000  1.00  0.00           N
ATOM      7  CA  CYS A  13       1.900  -2.000   0.000  1.00  0.00           C
ATOM      8  C   CYS A  13       3.000  -1.100   0.000  1.00  0.00           C
ATOM      9  O   CYS A  13       4.100  -1.500   0.000  1.00  0.00           O
ATOM     10  SG  CYS A  13       2.200   0.000   0.000  1.00  0.00           S
HETATM   11 ZN    ZN A 101       0.000   0.000   0.000  1.00  0.00          ZN
TER
END
"""
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "tiny.pdb"
        path.write_text(pdb)
        result = process_structure(
            path,
            ProcessingOptions(
                fill_missing_residues=False,
                add_missing_atoms=False,
                amber_first_phosphate_cleanup=False,
            ),
        )
    assert result.output_format == "pdb"
    assert "HIE A  12" in result.output_text
    assert "CYM A  13" in result.output_text
    print("smoke workflow ok")


if __name__ == "__main__":
    main()
