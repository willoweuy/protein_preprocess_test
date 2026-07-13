from pathlib import Path
import json
import subprocess
import sys
from tempfile import TemporaryDirectory


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
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
    with TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        input_path = tmp / "input.pdb"
        output_path = tmp / "output.pdb"
        report_path = tmp / "report.json"
        config_path = tmp / "config.json"
        input_path.write_text(pdb)
        config_path.write_text(
            json.dumps(
                {
                    "fill_missing_residues": False,
                    "add_missing_atoms": False,
                    "amber_first_phosphate_cleanup": False,
                    "output_format": "pdb",
                }
            )
        )
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "protein_process.cli",
                "--input",
                str(input_path),
                "--output",
                str(output_path),
                "--config",
                str(config_path),
                "--report",
                str(report_path),
            ],
            cwd=str(project_root),
            capture_output=True,
            text=True,
            check=False,
        )
        assert completed.returncode == 0, completed.stderr or completed.stdout
        report = json.loads(report_path.read_text())
        output = output_path.read_text()
        assert report["status"] == "ok"
        assert "HIE A  12" in output
        assert "CYM A  13" in output
        print("smoke cli ok")


if __name__ == "__main__":
    main()
