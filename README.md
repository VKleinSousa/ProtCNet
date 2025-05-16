# ProtCNet: Protein Contact Network Analyzer

**ProtCNet** is an interactive Python tool (Colab notebook) for analyzing and visualizing residue- and chain-level contact networks in protein structures. 

---

## 🔍 What It Does

* Parses **.pdb** or **.cif** files (upload or download via PDB ID)
* Selects atoms of interest (e.g., **Cα**, **Cβ**, or all atoms)
* Filters residues by interaction type (e.g., **hydrophobic**, **electrostatic**)
* Computes:

  * **Inter-chain** contacts
  * **Intra-chain** contacts with sequence distance control
  * **All** contacts (hybrid)
* Builds and visualizes:

  * Chain-level networks (**interactive HTML**)
  * Residue-level networks (**interactive HTML**)
  * Self-loop and detailed edge diagrams (**high-res TIFF**)
  * Heatmaps of contact frequency
* Exports **contact tables** for downstream use in **ChimeraX**, analysis, or annotation

---
## 🖼️ Example Visualizations
# Inter-chain contacts
![Inter-chain level contacts](https://github.com/VKleinSousa/ProtCNet/blob/main/5IV5_chainlevelnet.png)

# Inter-chain residue-level contacts
![Residue-level contacts](https://github.com/VKleinSousa/ProtCNet/blob/main/5IV5_residuelevelnet.png)

## 🚀 Getting Started

## 🚀 Getting Started

### Option 1: Run via Google Colab (Recommended)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1NCpjE7R8TzlSWjKSaQOj-lR_OvJuBUEc?usp=sharing)


### Option 2: Clone and run locally

```bash
git clone https://github.com/VKleinSousa/ProtCNet/
cd ProtCNet
jupyter notebook ProtCNet.ipynb
```

---

## 🌐 Output Files

* `network_visualization.html`: Chain-level network (interactive)
* `network_residue_level.html`: Residue-level network (interactive)
* `network_visualization.tiff`: High-resolution network with self-loops
* `heatmap.tiff`: Contact frequency heatmap
* `*_contacts.tsv`: Tab-delimited contact table (ChimeraX compatible)

---

## 🚚 Dependencies

* `biopython`
* `networkx`
* `scipy`
* `seaborn`
* `matplotlib`
* `ipywidgets`
* `ipysigma`

Install with:

```bash
pip install biopython networkx scipy seaborn matplotlib ipywidgets ipysigma
```

---


## 🌐 Contributors

Developed by Victor Klein-Sousa (@vkleinsousa) under the MIT License.

If you use this tool for your research, please cite the GitHub repository.

GitHub: [https://github.com/VKleinSousa/ProtCNet](https://github.com/VKleinSousa/ProtCNet)

---
