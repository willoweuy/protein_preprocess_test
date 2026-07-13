from __future__ import annotations

import json
from pathlib import Path
import tempfile

import streamlit as st
import streamlit.components.v1 as components

from protein_process import (
    LongGapDecisionRequired,
    MutationRequest,
    ProcessingOptions,
    process_structure,
)
from protein_process.pdb_text import apply_altloc_selections, find_altloc_residues


st.set_page_config(page_title="Protein Process", layout="wide")


def _write_uploaded_file(uploaded_file, suffix: str) -> Path:
    handle = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    handle.write(uploaded_file.getbuffer())
    handle.flush()
    handle.close()
    return Path(handle.name)


def _write_text_file(text: str, suffix: str) -> Path:
    handle = tempfile.NamedTemporaryFile("w+", delete=False, suffix=suffix)
    handle.write(text)
    handle.flush()
    handle.close()
    return Path(handle.name)


def _suffix_for_upload(uploaded_file, fallback: str) -> str:
    name = uploaded_file.name.lower()
    for suffix in (".cif", ".mmcif", ".pdb", ".fa", ".fasta"):
        if name.endswith(suffix):
            return suffix
    return fallback


def _uploaded_text(uploaded_file) -> str:
    return uploaded_file.getvalue().decode("utf-8", errors="replace")


def _mutation_rows_to_requests(rows) -> tuple[MutationRequest, ...]:
    requests: list[MutationRequest] = []
    for row in rows:
        chain_id = str(row.get("chain_id", "")).strip()
        residue_number = str(row.get("residue_number", "")).strip()
        insertion_code = str(row.get("insertion_code", "")).strip()
        new_resname = str(row.get("new_resname", "")).strip().upper()
        if not chain_id and not residue_number and not new_resname:
            continue
        if not chain_id or not residue_number or not new_resname:
            st.warning("Skipped one mutation row because chain_id, residue_number, and new_resname are all required.")
            continue
        requests.append(
            MutationRequest(
                chain_id=chain_id,
                residue_number=residue_number,
                insertion_code=insertion_code,
                new_resname=new_resname,
            )
        )
    return tuple(requests)


def _mutation_target_options(rows) -> list[tuple[str, str, str, str]]:
    options: list[tuple[str, str, str, str]] = []
    for row in rows:
        chain_id = str(row.get("chain_id", "")).strip()
        residue_number = str(row.get("residue_number", "")).strip()
        insertion_code = str(row.get("insertion_code", "")).strip()
        new_resname = str(row.get("new_resname", "")).strip().upper()
        if chain_id and residue_number:
            label = f"Mutation target: chain {chain_id} residue {residue_number}{insertion_code} {new_resname}"
            options.append((label, chain_id, residue_number, insertion_code))
    return options


def _render_structure_viewer(
    structure_text: str,
    structure_format: str,
    chain_id: str,
    residue_number: str,
    height: int = 560,
) -> None:
    if not structure_text.strip():
        st.info("No structure text available for the viewer yet.")
        return

    selection = {}
    if chain_id.strip() and residue_number.strip():
        selection["chain"] = chain_id.strip()
        try:
            selection["resi"] = int(residue_number)
        except ValueError:
            selection["resi"] = residue_number.strip()

    html = f"""
<div id="protein-viewer" style="height:{height}px; width:100%; position:relative;"></div>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script>
const data = {json.dumps(structure_text)};
const format = {json.dumps(structure_format)};
const target = {json.dumps(selection)};
const viewer = $3Dmol.createViewer("protein-viewer", {{backgroundColor: "white"}});
viewer.addModel(data, format);
viewer.setStyle({{}}, {{cartoon: {{color: "spectrum"}}, stick: {{radius: 0.12}}}});
if (Object.keys(target).length > 0) {{
  viewer.addStyle(target, {{stick: {{radius: 0.35, colorscheme: "greenCarbon"}}, sphere: {{radius: 0.45, color: "lime"}}}});
  viewer.zoomTo(target);
}} else {{
  viewer.zoomTo();
}}
viewer.render();
</script>
"""
    components.html(html, height=height + 20)


def _render_long_gap_decision() -> str:
    gaps = st.session_state.get("long_gap_gaps", [])
    if not gaps:
        return "ask"

    st.warning("One or more missing internal residue segments are longer than the threshold.")
    for gap in gaps:
        st.write(
            f"Chain `{gap.chain_id or '_'}` after residue `{gap.after_residue or '?'}`: "
            f"{gap.length} missing residues."
        )
    choice = st.radio(
        "How should the run proceed?",
        [
            "Stop until I decide",
            "Continue filling long missing segments",
            "Treat long gaps as split pieces and add SEGID labels",
        ],
        index=0,
    )
    if choice == "Continue filling long missing segments":
        return "continue"
    if choice == "Treat long gaps as split pieces and add SEGID labels":
        return "split"
    return "ask"


st.title("Interactive Protein Process")

structure_file = st.file_uploader("Protein structure", type=["cif", "mmcif", "pdb"])
fasta_file = st.file_uploader("Optional FASTA for PDB + FASTA workflows", type=["fa", "fasta"])

long_gap_action = _render_long_gap_decision()

uploaded_structure_text = _uploaded_text(structure_file) if structure_file is not None else ""
structure_suffix = _suffix_for_upload(structure_file, ".pdb") if structure_file is not None else ".pdb"
is_pdb_upload = structure_suffix == ".pdb"
altloc_residues = find_altloc_residues(uploaded_structure_text) if is_pdb_upload and uploaded_structure_text else ()

repair_tab, force_field_tab, occupancy_tab, mutate_tab, viewer_tab, download_tab = st.tabs(
    ["Repair", "Force Field", "Occupancy", "Mutate", "Viewer", "Download"]
)

with repair_tab:
    st.subheader("Structure Repair")
    fill_missing_residues = st.checkbox("Fill missing internal residues", value=True)
    add_missing_atoms = st.checkbox("Fill missing atoms", value=True)
    add_terminal_caps = st.checkbox("Add ACE/NME terminal caps", value=True)
    add_gromacs_ter_records = st.checkbox("Add GROMACS TER molecule separation", value=True)
    add_hydrogens = st.checkbox("Add missing hydrogens", value=False)
    ph = st.number_input("pH for hydrogens and PropKa decisions", min_value=0.0, max_value=14.0, value=7.4, step=0.1)
    long_gap_threshold = st.number_input(
        "Warn before filling a missing segment longer than this residue count",
        min_value=1,
        max_value=200,
        value=20,
        step=1,
    )

with force_field_tab:
    st.subheader("Zn and Amber Options")
    zn_cutoff = st.slider("Zn coordination cutoff (A)", min_value=1.8, max_value=3.5, value=2.8, step=0.1)
    histidine_mode = st.selectbox(
        "Zn-coordinating histidine naming",
        ["auto", "HID", "HIE", "HIP"],
        index=0,
        help="Auto uses ND1 -> HID, NE2 -> HIE, and both donors -> HIP.",
    )
    amber_cleanup = st.checkbox("Amber DNA/RNA first-residue P/O1P/O2P/OP1/OP2 cleanup", value=True)

with occupancy_tab:
    st.subheader("Alternate Occupancy")
    altloc_selections: dict[str, str] = {}
    if structure_file is None:
        st.info("Upload a PDB file to inspect alternate occupancy records.")
    elif not is_pdb_upload:
        st.info("Manual alternate-location selection is currently available for PDB uploads.")
    elif not altloc_residues:
        st.info("No alternate-location records were detected in this PDB upload.")
    else:
        st.write("Default selection is highest mean occupancy; exact ties choose `A` when present.")
        for residue in altloc_residues:
            option_labels = {
                choice: f"{choice} (mean occupancy {residue.occupancy_by_choice.get(choice, 0.0):.2f})"
                for choice in residue.choices
            }
            label = (
                f"Chain {residue.chain_id or '_'} residue "
                f"{residue.residue_number}{residue.insertion_code} {residue.resname}"
            )
            altloc_selections[residue.key] = st.selectbox(
                label,
                residue.choices,
                index=residue.choices.index(residue.default_choice),
                format_func=lambda choice, labels=option_labels: labels[choice],
                key=f"altloc_{residue.key}",
            )

with mutate_tab:
    st.subheader("Manual Mutation and Protonation")
    use_propka = st.checkbox("Use PropKa to predict protonation state", value=False)
    propka_apply = st.checkbox("Apply PropKa residue-name suggestions", value=False, disabled=not use_propka)
    st.caption("Add rows like chain A, residue 12, new_resname HIP. Protonation labels are applied as residue-name changes.")
    mutation_rows = st.data_editor(
        [{"chain_id": "", "residue_number": "", "insertion_code": "", "new_resname": ""}],
        num_rows="dynamic",
        use_container_width=True,
        column_config={
            "chain_id": st.column_config.TextColumn("chain_id", max_chars=2),
            "residue_number": st.column_config.TextColumn("residue_number"),
            "insertion_code": st.column_config.TextColumn("insertion_code", max_chars=1),
            "new_resname": st.column_config.TextColumn("new_resname", max_chars=3),
        },
    )

with viewer_tab:
    st.subheader("Structure Viewer")
    result = st.session_state.get("result")
    viewer_source_options = ["Uploaded structure"]
    if result is not None:
        viewer_source_options.insert(0, "Processed structure")
    viewer_source = st.radio("Viewer source", viewer_source_options, horizontal=True)

    target_options = [("Manual target", "", "", "")]
    target_options.extend(_mutation_target_options(mutation_rows))
    for residue in altloc_residues:
        target_options.append(
            (
                f"Occupancy target: chain {residue.chain_id or '_'} residue {residue.residue_number}{residue.insertion_code} {residue.resname}",
                residue.chain_id,
                residue.residue_number,
                residue.insertion_code,
            )
        )
    selected_target = st.selectbox("Zoom target", target_options, format_func=lambda item: item[0])

    if selected_target[0] == "Manual target":
        col1, col2, _col3 = st.columns([1, 1, 1])
        with col1:
            viewer_chain = st.text_input("Viewer chain", value="")
        with col2:
            viewer_residue = st.text_input("Viewer residue number", value="")
    else:
        viewer_chain = selected_target[1]
        viewer_residue = selected_target[2]
        st.write(f"Zooming to chain `{viewer_chain or '_'}` residue `{viewer_residue}`.")

    if viewer_source == "Processed structure" and result is not None:
        viewer_text = result.output_text
        viewer_format = "pdb" if result.output_format == "pdb" else "mmcif"
    else:
        viewer_text = uploaded_structure_text
        viewer_format = "pdb" if is_pdb_upload else "mmcif"
    _render_structure_viewer(viewer_text, viewer_format, viewer_chain, viewer_residue)

with download_tab:
    st.subheader("Download")
    output_format = st.radio("Download format", ["pdb", "cif"], horizontal=True)

run = st.button("Run", type="primary", disabled=structure_file is None)

if run and structure_file is not None:
    altloc_report: list[str] = []
    if is_pdb_upload:
        selected_structure_text, altloc_report = apply_altloc_selections(
            uploaded_structure_text,
            altloc_selections,
        )
        structure_path = _write_text_file(selected_structure_text, ".pdb")
    else:
        structure_path = _write_uploaded_file(structure_file, structure_suffix)
    fasta_path = None
    if fasta_file is not None:
        fasta_path = _write_uploaded_file(fasta_file, _suffix_for_upload(fasta_file, ".fasta"))

    options = ProcessingOptions(
        fill_missing_residues=fill_missing_residues,
        add_missing_atoms=add_missing_atoms,
        add_hydrogens=add_hydrogens,
        ph=float(ph),
        long_gap_threshold=int(long_gap_threshold),
        long_gap_action=long_gap_action,
        zn_cutoff_angstrom=float(zn_cutoff),
        histidine_mode=histidine_mode,
        amber_first_phosphate_cleanup=amber_cleanup,
        add_terminal_caps=add_terminal_caps,
        add_gromacs_ter_records=add_gromacs_ter_records,
        use_propka=use_propka,
        propka_apply_predictions=propka_apply,
        mutations=_mutation_rows_to_requests(mutation_rows),
        output_format=output_format,
    )

    try:
        result = process_structure(structure_path, options, fasta_path=fasta_path)
        result.report = altloc_report + result.report
    except LongGapDecisionRequired as exc:
        st.session_state.long_gap_gaps = exc.gaps
        st.error("Run paused. Choose how to handle the long missing segment warning, then click Run again.")
        st.stop()
    except RuntimeError as exc:
        st.error(str(exc))
        st.stop()
    except Exception as exc:
        st.exception(exc)
        st.stop()

    st.session_state.result = result
    st.session_state.long_gap_gaps = []

result = st.session_state.get("result")
if result is not None:
    st.success("Processing complete.")
    if result.warnings:
        st.subheader("Warnings")
        for warning in result.warnings:
            st.warning(warning)
    st.subheader("Report")
    st.code("\n".join(result.report), language="text")
    filename = f"processed.{result.output_format}"
    st.download_button(
        "Download processed structure",
        data=result.output_text,
        file_name=filename,
        mime="chemical/x-pdb" if result.output_format == "pdb" else "chemical/x-cif",
    )
