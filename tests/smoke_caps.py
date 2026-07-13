from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from protein_process.pdb_text import add_terminal_caps, get_segid, parse_atom_line, set_segid


def _atom(serial: int, atom_name: str, resseq: int, x: float, resname: str = "ALA") -> str:
    element = "".join(ch for ch in atom_name if ch.isalpha())[:1]
    return (
        f"ATOM  {serial:5d} {atom_name:>4} {resname:>3} A{resseq:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           {element:>2}\n"
    )


def main() -> None:
    pdb = "".join(
        [
            set_segid(_atom(1, "N", 10, 0.0), "D1"),
            set_segid(_atom(2, "CA", 10, 1.0), "D1"),
            set_segid(_atom(3, "C", 10, 2.0), "D1"),
            set_segid(_atom(4, "O", 10, 3.0), "D1"),
            set_segid(_atom(5, "OXT", 10, 3.2), "D1"),
        ]
    )
    capped, report, warnings = add_terminal_caps(pdb)
    assert not warnings
    assert "ACE" in capped
    assert "NME" in capped
    assert "OXT" not in capped
    assert "ALA A  10" in capped
    assert "ACE A   9" in capped
    assert "NME A  11" in capped
    cap_lines = [line for line in capped.splitlines(keepends=True) if " ACE " in line or " NME " in line]
    assert cap_lines
    assert all(get_segid(line) == "D1" for line in cap_lines), capped
    atoms = [parse_atom_line(i, line) for i, line in enumerate(capped.splitlines(keepends=True))]
    atoms = [atom for atom in atoms if atom is not None]
    by_res_atom = {(atom.resname, atom.atom_name): atom for atom in atoms}
    ace_c = by_res_atom[("ACE", "C")]
    ace_ch3 = by_res_atom[("ACE", "CH3")]
    ace_o = by_res_atom[("ACE", "O")]
    ala_n = by_res_atom[("ALA", "N")]
    ala_ca = by_res_atom[("ALA", "CA")]
    ala_c = by_res_atom[("ALA", "C")]
    nma_n = by_res_atom[("NME", "N")]
    nma_ch3 = by_res_atom[("NME", "CH3")]
    assert 1.25 < _distance(ace_c, ala_n) < 1.40
    assert 1.15 < _distance(ace_c, ace_o) < 1.30
    assert 1.40 < _distance(ace_c, ace_ch3) < 1.60
    assert 1.25 < _distance(ala_c, nma_n) < 1.40
    assert 1.35 < _distance(nma_n, nma_ch3) < 1.55
    assert 115.0 < _angle(ala_ca, ala_n, ace_c) < 125.0
    assert 115.0 < _angle(ala_c, nma_n, nma_ch3) < 125.0
    assert "added ACE" in "\n".join(report)
    assert "added NME" in "\n".join(report)
    print("smoke caps ok")


def _distance(atom_a, atom_b) -> float:
    return ((atom_a.x - atom_b.x) ** 2 + (atom_a.y - atom_b.y) ** 2 + (atom_a.z - atom_b.z) ** 2) ** 0.5


def _angle(atom_a, atom_b, atom_c) -> float:
    import math

    ba = (atom_a.x - atom_b.x, atom_a.y - atom_b.y, atom_a.z - atom_b.z)
    bc = (atom_c.x - atom_b.x, atom_c.y - atom_b.y, atom_c.z - atom_b.z)
    ba_length = (ba[0] ** 2 + ba[1] ** 2 + ba[2] ** 2) ** 0.5
    bc_length = (bc[0] ** 2 + bc[1] ** 2 + bc[2] ** 2) ** 0.5
    cosine = (ba[0] * bc[0] + ba[1] * bc[1] + ba[2] * bc[2]) / (ba_length * bc_length)
    cosine = min(1.0, max(-1.0, cosine))
    return math.degrees(math.acos(cosine))


if __name__ == "__main__":
    main()
