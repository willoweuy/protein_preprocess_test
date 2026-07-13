from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

from .models import MutationRequest


def run_propka_prediction(pdb_path: Path, ph: float) -> tuple[list[MutationRequest], list[str]]:
    executable = shutil.which("propka3") or shutil.which("propka")
    if executable is None:
        return [], ["PropKa was selected, but neither `propka3` nor `propka` was found on PATH."]

    try:
        completed = subprocess.run(
            [executable, str(pdb_path)],
            cwd=str(pdb_path.parent),
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except Exception as exc:
        return [], [f"PropKa failed to start: {exc}"]

    report = [f"PropKa command: {executable} {pdb_path.name}"]
    if completed.returncode != 0:
        stderr = completed.stderr.strip()
        if stderr:
            report.append(f"PropKa returned non-zero status: {stderr}")
        else:
            report.append("PropKa returned non-zero status with no stderr output.")
        return [], report

    pka_path = pdb_path.with_suffix(".pka")
    if not pka_path.exists():
        report.append("PropKa finished, but no .pka file was produced.")
        return [], report

    mutations: list[MutationRequest] = []
    summary_line = re.compile(
        r"^\s*(ASP|GLU|HIS|LYS|CYS)\s+(-?\d+)\s+([A-Za-z0-9])\s+(-?\d+(?:\.\d+)?)"
    )
    for line in pka_path.read_text(errors="replace").splitlines():
        match = summary_line.match(line)
        if not match:
            continue
        resname, resseq, chain_id, pka_value = match.groups()
        pka = float(pka_value)
        new_resname: str | None = None
        if resname == "ASP" and pka > ph:
            new_resname = "ASH"
        elif resname == "GLU" and pka > ph:
            new_resname = "GLH"
        elif resname == "HIS" and pka > ph:
            new_resname = "HIP"
        elif resname == "LYS" and pka < ph:
            new_resname = "LYN"
        elif resname == "CYS" and pka < ph:
            new_resname = "CYM"
        if new_resname:
            mutations.append(MutationRequest(chain_id=chain_id, residue_number=resseq, new_resname=new_resname))
            report.append(
                f"PropKa pKa {pka:.2f}: chain {chain_id} residue {resseq} {resname} -> suggested {new_resname} at pH {ph:.2f}."
            )

    if not mutations:
        report.append(f"PropKa produced no automatic residue-name suggestions at pH {ph:.2f}.")
    return mutations, report
