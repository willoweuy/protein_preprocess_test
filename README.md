# Interactive Protein Process

Local Streamlit module for first-pass protein structure preprocessing. This version supports CIF or PDB input, optional FASTA upload, missing residue and atom repair through PDBFixer/OpenMM, ACE/NME terminal capping, PDB alternate-location occupancy selection, Zn coordination residue naming, Amber first-nucleotide phosphate cleanup, manual residue-name changes, optional PropKa reporting, visualization, and PDB/CIF export.

Small ligand preparation is intentionally out of scope for this first version.

## Setup

```bash
conda env create -p ./.conda-env -f environment.yml
conda activate ./.conda-env
streamlit run app.py
```

If you prefer an existing environment, install the packages from `requirements.txt`. The conda route is recommended because `pdbfixer` and `openmm` are most reliable from conda-forge.

## Tests

Run the smoke tests from the repository root:

```bash
python tests/smoke_workflow.py
python tests/smoke_cli.py
python tests/smoke_fasta.py
python tests/smoke_altloc.py
python tests/smoke_amber_cleanup.py
python tests/smoke_caps.py
python tests/smoke_gromacs_ter.py
python tests/smoke_mmcif_ids.py
python tests/smoke_segid.py
```

## Workflow

1. Upload a `.cif`, `.mmcif`, or `.pdb` structure. Optionally upload FASTA for PDB plus FASTA workflows.
2. Select repair options, Zn naming options, Amber cleanup, manual mutation/protonation rows, PropKa options, and output format.
3. For PDB files with alternate-location records, review the **Occupancy** tab. Defaults choose the highest mean occupancy; exact ties choose `A` when present.
4. Use the **Viewer** tab to inspect the uploaded or processed structure and zoom to mutation or occupancy targets.
5. Click **Run**.
6. If PDBFixer reports a missing internal segment longer than the threshold, the run pauses. Choose whether to continue filling or treat the gap as split pieces with SEGID labels, then click **Run** again.
7. Download the processed structure as PDB or CIF.

## Notes

- Default pH is 7.4.
- Default long-gap threshold is 20 residues.
- ACE/NME terminal capping is enabled by default. Existing protein residue numbers are preserved; caps are added as adjacent residues where possible.
- GROMACS TER molecule separation is enabled by default. It inserts `TER` records between same-chain molecule-type boundaries, especially around ligands/ions/waters, without changing chain IDs, SEGIDs, or residue numbers.
- Terminal missing residues reported by PDBFixer are ignored before the long-gap warning and filling step.
- Default PropKa behavior is off.
- PDB alternate-location selection is applied before PDBFixer sees the structure. Manual selection currently appears for PDB uploads with explicit altLoc records.
- Default Amber cleanup removes `P`, `O1P`, `O2P`, `OP1`, and `OP2` from the smallest-numbered DNA/RNA residue in each nucleic-acid chain. If more than 3 of those atoms are removed from one first residue, the workflow emits a warning.
- Normal output preserves existing SEGID fields and chain IDs.
- Split-gap output uses numbered SEGIDs from the chain ID, such as `D1`, `D2`, to avoid collisions with real neighboring chains.
- For CIF/mmCIF inputs, output PDB chain IDs are restored from `_atom_site.auth_asym_id`, and PDB SEGID is restored from `_atom_site.label_asym_id` when available. This keeps PyMOL-style identifiers closer to the original loaded mmCIF.
- Zn coordination naming converts coordinating `CYS` to `CYM`; coordinating histidines are set to `HID`, `HIE`, or `HIP` based on donor atoms unless you choose a fixed histidine name.
- Manual protonation labels such as `HID`, `HIE`, `HIP`, `ASH`, `GLH`, `CYM`, and `LYN` are applied as residue-name changes.
- Standard amino-acid mutations are attempted through PDBFixer first. If PDBFixer cannot rebuild the mutation, the app falls back to residue-name replacement and reports a warning.

## Known Limitations

- Small-ligand preparation and ligand parameterization are intentionally out of scope.
- Missing-residue filling depends on sequence metadata. Original mmCIF/PDB files or matching FASTA files are more reliable than coordinates exported from PyMOL.
- ACE/NME cap placement uses a local heuristic with trigonal geometry checks. It is intended to avoid obviously bad geometry, not to replace full minimization.
- PropKa support depends on the local `propka` installation and is off by default.
- Manual alternate-location selection is currently focused on PDB altLoc records.

## Repository Notes

- `references/add_caps_reference.py` is preserved as reference material from development. The active capping implementation lives in `protein_process/pdb_text.py`.
- Local environments, caches, generated PropKa files, and macOS metadata are ignored by `.gitignore`.

## PyMOL Plugin

A PyMOL plugin is available at:

```text
pymol_plugin/protein_preprocess_plugin.py
```

Install it from **PyMOL > Plugin > Plugin Manager > Install New Plugin**. The plugin can process the selected PyMOL object export with an optional sequence file (`.fa`, `.fasta`, `.pdb`, `.cif`, `.mmcif`) for missing internal residues. It remembers the project/Python paths, can download the processed PDB or CIF to a selected output path with overwrite confirmation, and can save the full report as a log file. See `pymol_plugin/README.md` for details.
