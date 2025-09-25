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
* Exports **contact tables** for downstream use.

---
## 🖼️ Example Visualizations
# Inter-chain contacts
![Inter-chain level contacts](https://github.com/VKleinSousa/ProtCNet/blob/main/5IV5_chainlevelnet.png)

# Inter-chain residue-level contacts
![Residue-level contacts](https://github.com/VKleinSousa/ProtCNet/blob/main/5IV5_residuelevelnet.png)


## 🚀 Getting Started

### Option 1: Run via Google Colab (Recommended)

Version 0.2:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1_-mGMYDDJY6na9hQP21KO_4Bwhi8tXEO?usp=sharing)

Version 0.1

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1NCpjE7R8TzlSWjKSaQOj-lR_OvJuBUEc?usp=sharing)


### Option 2: Clone and run locally

```bash
git clone https://github.com/VKleinSousa/ProtCNet/
cd ProtCNet
jupyter notebook ProtCNet.ipynb
```

**Recommended:**
Running protcnet_cli.py locally:

```bash
conda create -n protcnet \
  python=3.10 \
  biopython \
  networkx \
  matplotlib \
  seaborn \
  scipy \
  pip \
  -c conda-forge -c bioconda \
  mmseqs2 -y

conda activate protcnet

pip install ipysigma
```

Finally, you can run:
```
python protcnet_cli.py -h
```

```
usage: protcnet_cli.py [-h] [--pdb_id PDB_ID] [--input_file INPUT_FILE]
                       [--exclude_fasta EXCLUDE_FASTA]
                       [--exclude_chain EXCLUDE_CHAIN [EXCLUDE_CHAIN ...]]
                       [--cutoff CUTOFF] [--mode {inter,intra,all}]
                       [--sequence_distance SEQUENCE_DISTANCE]
                       [--atom_mode {all,ca,cb}]
                       [--interaction_filter {all,electrostatic,hydrophobic}]
                       [--track_residues] [--cluster_chains]
                       [--weight_contacts] [--residue_level_net] [--force]
                       [--residue_level_cutoff RESIDUE_LEVEL_CUTOFF]
                       [--network_file NETWORK_FILE]
                       [--edge_style {rectangle,line,curve}]

ProtCNet CLI (refactored from Colab)

options:
  -h, --help            show this help message and exit
  --pdb_id PDB_ID       PDB ID to download
  --input_file INPUT_FILE
                        Path to local .pdb or .cif file
  --exclude_fasta EXCLUDE_FASTA
                        FASTA file with sequences to exclude if similar (>0.8)
  --exclude_chain EXCLUDE_CHAIN [EXCLUDE_CHAIN ...]
                        List of chain IDs to exclude from contact network
                        (e.g. A B C)
  --cutoff CUTOFF       Distance cutoff (Å)
  --mode {inter,intra,all}
                        Contact type
  --sequence_distance SEQUENCE_DISTANCE
                        Minimum residue separation for intra-chain contacts
  --atom_mode {all,ca,cb}
                        Atoms to use for contact calculation
  --interaction_filter {all,electrostatic,hydrophobic}
                        Residue filtering
  --track_residues      Track contacting residues
  --cluster_chains      Cluster chains by sequence similarity (MMseqs2)
  --weight_contacts     Scale node size by contact count
  --residue_level_net   Plot residue-level contact network
  --force               Force overun mmseqs2, will delete previous mmseqs2
                        run.
  --residue_level_cutoff RESIDUE_LEVEL_CUTOFF
                        Minimum contacts for residue-level node display
  --network_file NETWORK_FILE
                        Output file for network
  --edge_style {rectangle,line,curve}
                        Edge style in network
```

---
Example:

``` python protcnet_cli.py --input_file pdbs/5IV5.cif --cutoff 8.0 --cluster_chains --edge_style rectangle --network_file T4.html --force --exclude_fasta exclude_fasta.fa --exclude_chain BB BC DE DF FH FI IA IB R S o p BG DJ GC IF W t YA YB YC ZA ```


[🔗 View T4.html interactive network](https://vkleinsousa.github.io/ProtCNet/T4.html)


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

[![DOI](https://zenodo.org/badge/984670633.svg)](https://doi.org/10.5281/zenodo.15440764)

GitHub: [https://github.com/VKleinSousa/ProtCNet](https://github.com/VKleinSousa/ProtCNet)

---
