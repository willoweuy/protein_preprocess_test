# PyMOL Plugin

This plugin adds a **Plugin > Protein Preprocess** panel in PyMOL. It exports the selected PyMOL object as mmCIF, runs the project processor through the local conda environment, then loads the processed PDB or CIF back into PyMOL.

## Install

1. Make sure the project conda environment exists:

   ```bash
   conda env create -p ./.conda-env -f environment.yml
   ```

2. In PyMOL, open **Plugin > Plugin Manager > Install New Plugin**.
3. Select:

   ```text
   /path/to/protein_preprocess_test/pymol_plugin/protein_preprocess_plugin.py
   ```

4. Open **Plugin > Protein Preprocess**.

## Use

1. Load a structure into PyMOL.
2. Open the plugin panel.
3. Choose the PyMOL object.
4. Set **Project folder** to:

   ```text
   /path/to/protein_preprocess_test
   ```

5. Set **Processor Python** to:

   ```text
   /path/to/protein_preprocess_test/.conda-env/bin/python
   ```

6. For internal missing-residue filling, choose a `.fa`, `.fasta`, `.pdb`, `.cif`, or `.mmcif` file in **Sequence file**. FASTA record IDs should match chain IDs when there is more than one chain. Original `.cif/.mmcif/.pdb` files remain the most faithful source because they include both chain metadata and sequence metadata.
7. Leave **Remember paths** checked so these paths are restored next time.
8. Optionally set **Download file** and **Download format** to write the processed PDB or CIF to a permanent file. If blank, the result is loaded into PyMOL but the temporary file is not kept. If the selected file already exists, the plugin asks before overwriting it.
9. Choose repair, ACE/NME capping, GROMACS TER separation, Zn, Amber, PropKa, and mutation options.
10. Click **Run and Load**.

If a missing internal segment is longer than the threshold and the action is `ask`, the plugin prompts whether to continue filling or treat it as split pieces.

After completion, a resizable report window opens. Use **Save Log...** to write the full report to a `.log` or `.txt` file.

## Mutation Row Format

Use one mutation/protonation entry per line:

```text
A 12 HIP
A 12 B HIP
```

The four-column form includes insertion code.

## Notes

- The plugin loads the selected processed output format back into PyMOL.
- ACE/NME terminal capping is enabled by default. Existing protein residue numbers are preserved, cap residues are placed adjacent to the first and last protein residue numbers when possible, and cap bond placement uses 120-degree trigonal geometry at ACE `CA-N-C` and NME `C-N-CH3`.
- GROMACS TER separation is enabled by default. It inserts `TER` records between same-chain molecule-type boundaries so ligands and ions are not treated as one continuous chain by `pdb2gmx`.
- Alternate occupancy selection is applied automatically using the same default rule as the Streamlit app: highest mean occupancy, tie -> `A`.
- Detailed reports are shown after the processed object is loaded.
- If **Download file** is set, the plugin writes the selected processed PDB/CIF format there and records the path in the report.
- If **Sequence file** is blank, the plugin exports the selected PyMOL object to mmCIF and processes that. This is better than PDB as a fallback, but it may still lack the original polymer sequence metadata needed for internal missing-residue filling.
- If **Sequence file** is FASTA, the plugin exports the selected PyMOL object to mmCIF for coordinates and passes FASTA as PDBFixer's sequence source.
