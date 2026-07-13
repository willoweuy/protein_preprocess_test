from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import traceback

from pymol import cmd
from pymol.plugins import addmenuitemqt
from pymol.Qt import QtCore, QtWidgets


def _detect_project_root() -> Path:
    env_value = os.environ.get("PROTEIN_PREPROCESS_PROJECT", "")
    if env_value:
        env_path = Path(env_value).expanduser()
        if (env_path / "protein_process").exists():
            return env_path
    local_root = Path(__file__).resolve().parents[1]
    if (local_root / "protein_process").exists():
        return local_root
    return local_root


PLUGIN_ROOT = _detect_project_root()


def _default_python_path(project_root: Path) -> str:
    env_python = project_root / ".conda-env" / "bin" / "python"
    if env_python.exists():
        return str(env_python)
    return "python"


def _parse_mutations(text: str) -> list[dict]:
    mutations: list[dict] = []
    for line in text.splitlines():
        cleaned = line.strip()
        if not cleaned or cleaned.startswith("#"):
            continue
        parts = cleaned.replace(",", " ").split()
        if len(parts) == 3:
            chain_id, residue_number, new_resname = parts
            insertion_code = ""
        elif len(parts) == 4:
            chain_id, residue_number, insertion_code, new_resname = parts
        else:
            raise ValueError(
                f"Could not parse mutation line {line!r}. Use: A 12 HIP or A 12 B HIP."
            )
        mutations.append(
            {
                "chain_id": chain_id,
                "residue_number": residue_number,
                "insertion_code": insertion_code,
                "new_resname": new_resname.upper(),
            }
        )
    return mutations


def _format_report(payload: dict) -> str:
    lines: list[str] = []
    if payload.get("warnings"):
        lines.append("Warnings:")
        lines.extend(f"- {warning}" for warning in payload["warnings"])
        lines.append("")
    if payload.get("report"):
        lines.append("Report:")
        lines.extend(f"- {item}" for item in payload["report"])
    if not lines:
        lines.append(json.dumps(payload, indent=2))
    return "\n".join(lines)


class ReportDialog(QtWidgets.QDialog):
    def __init__(self, title: str, text: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(820, 620)
        layout = QtWidgets.QVBoxLayout(self)
        self.report_text = QtWidgets.QPlainTextEdit()
        self.report_text.setReadOnly(True)
        self.report_text.setPlainText(text)
        layout.addWidget(self.report_text, 1)

        buttons = QtWidgets.QDialogButtonBox()
        save_button = buttons.addButton("Save Log...", QtWidgets.QDialogButtonBox.ActionRole)
        close_button = buttons.addButton(QtWidgets.QDialogButtonBox.Close)
        save_button.clicked.connect(self._save_log)
        close_button.clicked.connect(self.accept)
        layout.addWidget(buttons)

    def _save_log(self):
        path, _selected_filter = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save Protein Preprocess Log",
            str(Path.home() / "protein_preprocess.log"),
            "Log files (*.log);;Text files (*.txt);;All files (*)",
        )
        if path:
            Path(path).write_text(self.report_text.toPlainText())


class ProteinPreprocessDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Protein Preprocess")
        self.resize(680, 760)
        self.settings = QtCore.QSettings("Codex", "ProteinPreprocess")

        layout = QtWidgets.QVBoxLayout(self)

        form = QtWidgets.QFormLayout()
        self.object_combo = QtWidgets.QComboBox()
        self._refresh_objects()
        refresh_button = QtWidgets.QPushButton("Refresh")
        refresh_button.clicked.connect(self._refresh_objects)
        object_row = QtWidgets.QHBoxLayout()
        object_row.addWidget(self.object_combo, 1)
        object_row.addWidget(refresh_button)
        form.addRow("PyMOL object", object_row)

        default_project_root = self.settings.value("project_root", str(PLUGIN_ROOT), type=str)
        if not (Path(default_project_root).expanduser() / "protein_process").exists():
            default_project_root = str(PLUGIN_ROOT)
        self.project_root = QtWidgets.QLineEdit(default_project_root)
        project_browse = QtWidgets.QPushButton("Browse")
        project_browse.clicked.connect(self._browse_project_root)
        project_row = QtWidgets.QHBoxLayout()
        project_row.addWidget(self.project_root, 1)
        project_row.addWidget(project_browse)
        form.addRow("Project folder", project_row)

        default_python = self.settings.value(
            "python_path",
            _default_python_path(Path(default_project_root)),
            type=str,
        )
        if default_python == "python" or not Path(default_python).expanduser().exists():
            default_python = _default_python_path(Path(default_project_root))
        self.python_path = QtWidgets.QLineEdit(default_python)
        python_browse = QtWidgets.QPushButton("Browse")
        python_browse.clicked.connect(self._browse_python_path)
        python_row = QtWidgets.QHBoxLayout()
        python_row.addWidget(self.python_path, 1)
        python_row.addWidget(python_browse)
        form.addRow("Processor Python", python_row)

        source_file_default = self.settings.value("source_file", "", type=str)
        self.source_file = QtWidgets.QLineEdit(source_file_default)
        source_browse = QtWidgets.QPushButton("Browse")
        source_browse.clicked.connect(self._browse_source_file)
        source_row = QtWidgets.QHBoxLayout()
        source_row.addWidget(self.source_file, 1)
        source_row.addWidget(source_browse)
        form.addRow("Sequence file", source_row)

        output_file_default = self.settings.value("output_file", "", type=str)
        self.output_file = QtWidgets.QLineEdit(output_file_default)
        output_browse = QtWidgets.QPushButton("Download / Save As...")
        output_browse.clicked.connect(self._browse_output_file)
        output_row = QtWidgets.QHBoxLayout()
        output_row.addWidget(self.output_file, 1)
        output_row.addWidget(output_browse)
        form.addRow("Download file", output_row)

        output_format_default = self.settings.value("output_format", "pdb", type=str).lower()
        self.output_format = QtWidgets.QComboBox()
        self.output_format.addItems(["pdb", "cif"])
        if output_format_default in {"pdb", "cif"}:
            self.output_format.setCurrentText(output_format_default)
        form.addRow("Download format", self.output_format)

        self.remember_paths = QtWidgets.QCheckBox("Remember paths")
        self.remember_paths.setChecked(self.settings.value("remember_paths", True, type=bool))
        form.addRow("", self.remember_paths)
        layout.addLayout(form)

        options_group = QtWidgets.QGroupBox("Options")
        options = QtWidgets.QFormLayout(options_group)
        self.fill_missing_residues = QtWidgets.QCheckBox()
        self.fill_missing_residues.setChecked(True)
        options.addRow("Fill internal missing residues", self.fill_missing_residues)

        self.add_missing_atoms = QtWidgets.QCheckBox()
        self.add_missing_atoms.setChecked(True)
        options.addRow("Fill missing atoms", self.add_missing_atoms)

        self.add_terminal_caps = QtWidgets.QCheckBox()
        self.add_terminal_caps.setChecked(True)
        options.addRow("Add ACE/NME terminal caps", self.add_terminal_caps)

        self.add_gromacs_ter_records = QtWidgets.QCheckBox()
        self.add_gromacs_ter_records.setChecked(True)
        options.addRow("Add GROMACS TER separation", self.add_gromacs_ter_records)

        self.add_hydrogens = QtWidgets.QCheckBox()
        options.addRow("Add hydrogens", self.add_hydrogens)

        self.ph = QtWidgets.QDoubleSpinBox()
        self.ph.setRange(0.0, 14.0)
        self.ph.setDecimals(1)
        self.ph.setSingleStep(0.1)
        self.ph.setValue(7.4)
        options.addRow("pH", self.ph)

        self.long_gap_threshold = QtWidgets.QSpinBox()
        self.long_gap_threshold.setRange(1, 200)
        self.long_gap_threshold.setValue(20)
        options.addRow("Long gap threshold", self.long_gap_threshold)

        self.long_gap_action = QtWidgets.QComboBox()
        self.long_gap_action.addItems(["ask", "continue", "split"])
        options.addRow("Long gap action", self.long_gap_action)

        self.histidine_mode = QtWidgets.QComboBox()
        self.histidine_mode.addItems(["auto", "HID", "HIE", "HIP"])
        options.addRow("Zn histidine naming", self.histidine_mode)

        self.zn_cutoff = QtWidgets.QDoubleSpinBox()
        self.zn_cutoff.setRange(1.8, 3.5)
        self.zn_cutoff.setDecimals(1)
        self.zn_cutoff.setSingleStep(0.1)
        self.zn_cutoff.setValue(2.8)
        options.addRow("Zn cutoff", self.zn_cutoff)

        self.amber_cleanup = QtWidgets.QCheckBox()
        self.amber_cleanup.setChecked(True)
        options.addRow("Amber first phosphate cleanup", self.amber_cleanup)

        self.use_propka = QtWidgets.QCheckBox()
        options.addRow("Use PropKa", self.use_propka)

        self.propka_apply = QtWidgets.QCheckBox()
        options.addRow("Apply PropKa suggestions", self.propka_apply)
        layout.addWidget(options_group)

        self.mutations = QtWidgets.QPlainTextEdit()
        self.mutations.setPlaceholderText("One per line: A 12 HIP\nOptional insertion code: A 12 B HIP")
        layout.addWidget(QtWidgets.QLabel("Manual mutation/protonation rows"))
        layout.addWidget(self.mutations, 1)

        buttons = QtWidgets.QDialogButtonBox()
        self.run_button = buttons.addButton("Run and Load", QtWidgets.QDialogButtonBox.AcceptRole)
        buttons.addButton(QtWidgets.QDialogButtonBox.Close)
        self.run_button.clicked.connect(self._run)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _browse_project_root(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select Protein Preprocess Project Folder",
            self.project_root.text(),
        )
        if path:
            self.project_root.setText(path)
            candidate_python = _default_python_path(Path(path))
            if candidate_python != "python":
                self.python_path.setText(candidate_python)

    def _browse_python_path(self):
        path, _selected_filter = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select Processor Python",
            self.python_path.text() or str(Path.home()),
            "Python executable (python python3);;All files (*)",
        )
        if path:
            self.python_path.setText(path)

    def _browse_source_file(self):
        path, _selected_filter = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select Sequence or Structure File",
            self.source_file.text() or str(Path.home()),
            "Sequence/structure files (*.fa *.fasta *.pdb *.cif *.mmcif);;All files (*)",
        )
        if path:
            self.source_file.setText(path)

    def _browse_output_file(self):
        output_format = self._selected_output_format()
        suffix = self._output_suffix()
        filters = "PDB files (*.pdb);;All files (*)" if output_format == "pdb" else "CIF files (*.cif);;All files (*)"
        path, _selected_filter = QtWidgets.QFileDialog.getSaveFileName(
            self,
            f"Download Processed {output_format.upper()} As",
            self.output_file.text() or str(Path.home() / f"processed{suffix}"),
            filters,
        )
        if path:
            if Path(path).suffix == "":
                path = f"{path}{suffix}"
            self.output_file.setText(path)

    def _selected_output_format(self) -> str:
        output_format = self.output_format.currentText().strip().lower()
        return output_format if output_format in {"pdb", "cif"} else "pdb"

    def _output_suffix(self) -> str:
        return ".cif" if self._selected_output_format() == "cif" else ".pdb"

    def _final_output_path(self) -> Path | None:
        output_file_text = self.output_file.text().strip()
        if not output_file_text:
            return None
        output_path = Path(output_file_text).expanduser().resolve()
        suffix = self._output_suffix()
        if output_path.suffix == "":
            output_path = output_path.with_suffix(suffix)
        elif output_path.suffix.lower() != suffix:
            output_path = output_path.with_suffix(suffix)
        return output_path

    def _confirm_overwrite(self, path: Path) -> bool:
        if not path.exists():
            return True
        answer = QtWidgets.QMessageBox.question(
            self,
            "Overwrite Existing File?",
            f"{path} already exists.\n\nOverwrite it?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.No,
        )
        return answer == QtWidgets.QMessageBox.Yes

    def _save_settings(self):
        self.settings.setValue("remember_paths", self.remember_paths.isChecked())
        self.settings.setValue("output_format", self._selected_output_format())
        if self.remember_paths.isChecked():
            self.settings.setValue("project_root", self.project_root.text())
            self.settings.setValue("python_path", self.python_path.text())
            self.settings.setValue("source_file", self.source_file.text())
            self.settings.setValue("output_file", self.output_file.text())

    def reject(self):
        self._save_settings()
        super().reject()

    def _refresh_objects(self):
        current = self.object_combo.currentText() if hasattr(self, "object_combo") else ""
        self.object_combo.clear()
        names = cmd.get_names("objects") or []
        self.object_combo.addItems(names)
        if current in names:
            self.object_combo.setCurrentText(current)

    def _config(self, long_gap_action: str | None = None) -> dict:
        return {
            "fill_missing_residues": self.fill_missing_residues.isChecked(),
            "add_missing_atoms": self.add_missing_atoms.isChecked(),
            "add_hydrogens": self.add_hydrogens.isChecked(),
            "ph": self.ph.value(),
            "long_gap_threshold": self.long_gap_threshold.value(),
            "long_gap_action": long_gap_action or self.long_gap_action.currentText(),
            "zn_cutoff_angstrom": self.zn_cutoff.value(),
            "histidine_mode": self.histidine_mode.currentText(),
            "amber_first_phosphate_cleanup": self.amber_cleanup.isChecked(),
            "add_terminal_caps": self.add_terminal_caps.isChecked(),
            "add_gromacs_ter_records": self.add_gromacs_ter_records.isChecked(),
            "use_propka": self.use_propka.isChecked(),
            "propka_apply_predictions": self.propka_apply.isChecked(),
            "mutations": _parse_mutations(self.mutations.toPlainText()),
            "output_format": self._selected_output_format(),
            "apply_altloc_selection": True,
        }

    def _run_cli(
        self,
        config: dict,
        input_path: Path,
        output_path: Path,
        report_path: Path,
        fasta_path: Path | None = None,
    ) -> tuple[int, dict]:
        config_path = report_path.with_suffix(".config.json")
        config_path.write_text(json.dumps(config, indent=2))
        project_root = Path(self.project_root.text()).expanduser().resolve()
        command = [
            self.python_path.text().strip(),
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
        ]
        if fasta_path is not None:
            command.extend(["--fasta", str(fasta_path)])
        completed = subprocess.run(
            command,
            cwd=str(project_root),
            capture_output=True,
            text=True,
            timeout=600,
        )
        payload = {}
        if report_path.exists():
            payload = json.loads(report_path.read_text())
        else:
            payload = {
                "status": "error",
                "message": completed.stderr.strip() or completed.stdout.strip() or "No report was produced.",
            }
        return completed.returncode, payload

    def _ask_long_gap_action(self, payload: dict) -> str | None:
        gap_lines = []
        for gap in payload.get("gaps", []):
            gap_lines.append(
                f"Chain {gap.get('chain_id') or '_'} after residue {gap.get('after_residue') or '?'}: "
                f"{gap.get('length')} missing residues"
            )
        box = QtWidgets.QMessageBox(self)
        box.setWindowTitle("Long Missing Segment")
        box.setText("One or more internal missing segments exceed the threshold.")
        box.setInformativeText("\n".join(gap_lines))
        continue_button = box.addButton("Continue filling", QtWidgets.QMessageBox.AcceptRole)
        split_button = box.addButton("Treat as split pieces", QtWidgets.QMessageBox.DestructiveRole)
        cancel_button = box.addButton(QtWidgets.QMessageBox.Cancel)
        box.exec_()
        clicked = box.clickedButton()
        if clicked == continue_button:
            return "continue"
        if clicked == split_button:
            return "split"
        if clicked == cancel_button:
            return None
        return None

    def _run(self):
        object_name = self.object_combo.currentText().strip()
        if not object_name:
            QtWidgets.QMessageBox.warning(self, "Protein Preprocess", "No PyMOL object selected.")
            return

        try:
            self._save_settings()
            with tempfile.TemporaryDirectory(prefix="protein_preprocess_") as tmpdir:
                tmp = Path(tmpdir)
                source_text = self.source_file.text().strip()
                sequence_path = Path(source_text).expanduser().resolve() if source_text else None
                fasta_path = None
                using_structure_sequence_file = False
                using_fasta_sequence_file = False
                if sequence_path is not None:
                    if not sequence_path.exists():
                        QtWidgets.QMessageBox.warning(
                            self,
                            "Protein Preprocess",
                            f"Sequence file does not exist:\n{sequence_path}",
                        )
                        return
                    suffix = sequence_path.suffix.lower()
                    if suffix in {".fa", ".fasta"}:
                        input_path = tmp / f"{object_name}.cif"
                        cmd.save(str(input_path), object_name, format="cif")
                        fasta_path = sequence_path
                        using_fasta_sequence_file = True
                    else:
                        input_path = sequence_path
                        using_structure_sequence_file = True
                else:
                    input_path = tmp / f"{object_name}.cif"
                    cmd.save(str(input_path), object_name, format="cif")
                output_format = self._selected_output_format()
                output_path = tmp / f"{object_name}_processed.{output_format}"
                report_path = tmp / "report.json"

                config = self._config()
                returncode, payload = self._run_cli(config, input_path, output_path, report_path, fasta_path)
                if returncode == 2 and payload.get("status") == "long_gap_required":
                    action = self._ask_long_gap_action(payload)
                    if action is None:
                        return
                    config = self._config(long_gap_action=action)
                    returncode, payload = self._run_cli(config, input_path, output_path, report_path, fasta_path)

                if returncode != 0 or payload.get("status") != "ok":
                    QtWidgets.QMessageBox.critical(
                        self,
                        "Protein Preprocess Failed",
                        payload.get("message", json.dumps(payload, indent=2)),
                    )
                    return

                if not using_structure_sequence_file and not using_fasta_sequence_file and self.fill_missing_residues.isChecked():
                    payload.setdefault("warnings", []).insert(
                        0,
                        "The plugin processed an mmCIF exported from PyMOL. This is a better fallback than PDB, "
                        "but PyMOL may still omit the original polymer sequence metadata needed to detect and fill "
                        "internal missing residues. Select a FASTA or original PDB/mmCIF sequence file in the plugin "
                        "for more reliable residue filling.",
                    )
                elif using_fasta_sequence_file and self.fill_missing_residues.isChecked():
                    payload.setdefault("warnings", []).insert(
                        0,
                        "The plugin used a FASTA sequence file with coordinates exported from PyMOL. This can detect "
                        "sequence gaps when FASTA record IDs match chain IDs, but an original mmCIF/PDB file remains "
                        "the most faithful source for chain metadata.",
                    )

                loaded_name = f"{object_name}_processed"
                final_output_path = self._final_output_path()
                if final_output_path is not None:
                    if self._confirm_overwrite(final_output_path):
                        final_output_path.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copyfile(output_path, final_output_path)
                        self.output_file.setText(str(final_output_path))
                        payload.setdefault("report", []).insert(
                            0,
                            f"Downloaded processed {output_format.upper()} to {final_output_path}.",
                        )
                    else:
                        payload.setdefault("warnings", []).insert(
                            0,
                            f"Skipped download because {final_output_path} already exists.",
                        )
                cmd.load(str(output_path), loaded_name)
                cmd.enable(loaded_name)
                cmd.zoom(loaded_name)
                ReportDialog("Protein Preprocess Complete", _format_report(payload), self).exec_()
        except Exception as exc:
            QtWidgets.QMessageBox.critical(
                self,
                "Protein Preprocess Error",
                f"{exc}\n\n{traceback.format_exc()}",
            )


def run_plugin_gui():
    dialog = ProteinPreprocessDialog()
    dialog.exec_()


def __init_plugin__(app=None):
    addmenuitemqt("Protein Preprocess", run_plugin_gui)
    cmd.extend("protein_preprocess_panel", run_plugin_gui)
