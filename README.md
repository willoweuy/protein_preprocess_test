# Protein Preprocess

Protein Preprocess is a local Streamlit app for preparing protein structures before molecular simulation. It accepts PDB or mmCIF input, can use an optional FASTA sequence file, repairs missing atoms or internal residues with PDBFixer/OpenMM, handles common zinc-coordination residue naming, applies simple protein terminal capping, and exports the processed structure as PDB or CIF.

Small-ligand preparation is not included in this version.

## Install

Create the recommended conda environment from the repository root:

```bash
conda env create -p ./.conda-env -f environment.yml
conda activate ./.conda-env
```

## Run the App

Start the Streamlit interface:

```bash
streamlit run app.py
```

Then open the local URL shown by Streamlit, upload a structure, choose the processing options, click **Run**, and download the processed PDB or CIF.

## PyMOL Plugin

A PyMOL plugin is included at:

```text
pymol_plugin/protein_preprocess_plugin.py
```

Install it from **PyMOL > Plugin > Plugin Manager > Install New Plugin**. In the plugin panel, select the PyMOL object, set the project folder and processor Python path, choose the processing options, then click **Run and Load**.

See [pymol_plugin/README.md](pymol_plugin/README.md) for plugin-specific setup and usage.
