#!/usr/bin/env python3
# ProtCNet CLI: Protein Contact Network Analyzer (Command-line version)

import os
import time
import argparse
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from scipy.spatial import KDTree
import csv
import subprocess
from ipysigma import Sigma
import shutil 

# ---------- Utility Functions ----------

def parse_structure(file_path, atom_mode="ca", residue_filter="all", with_residues=False):
    """
    Parse a PDB or CIF file and extract atomic coordinates grouped by chains.

    Returns:
        - chain_atoms: {chain_id: [atom_coords]}
        - atom_residue_map (optional): {chain_id: [(res_id, resname, resseq)]}
          where res_id = "A:101", resname = "GLY", resseq = 101
    """
    if file_path.endswith(".pdb"):
        parser = PDB.PDBParser(QUIET=True)
    elif file_path.endswith(".cif"):
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Only .pdb or .cif allowed.")

    structure = parser.get_structure("structure", file_path)

    chain_atoms = {}
    atom_residue_map = {}

    hydrophobic = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO'}
    charged = {'ASP', 'GLU', 'ARG', 'LYS', 'HIS'}

    for model in structure:
        for chain in model:
            chain_id = chain.id
            chain_atoms.setdefault(chain_id, [])
            if with_residues:
                atom_residue_map.setdefault(chain_id, [])

            for residue in chain:
                resname = residue.get_resname().strip()
                resseq = residue.id[1]  # sequence number
                res_id = f"{chain_id}:{resseq}"

                # Residue filtering
                if residue_filter == "hydrophobic" and resname not in hydrophobic:
                    continue
                if residue_filter == "electrostatic" and resname not in charged:
                    continue

                for atom in residue:
                    atom_name = atom.get_id()
                    if atom_mode == "ca" and atom_name != "CA":
                        continue
                    elif atom_mode == "cb" and atom_name != "CB":
                        continue

                    chain_atoms[chain_id].append(atom.coord)
                    if with_residues:
                        atom_residue_map[chain_id].append((res_id, resname, resseq))

    if with_residues:
        return chain_atoms, atom_residue_map
    else:
        return chain_atoms


def compute_contacts_with_residues(
    chain_atoms_dict,
    atom_residue_map,
    cutoff=5.0,
    mode="inter",
    sequence_distance=5,
    progress=None,
    status_label=None
):
    """
    Compute contacts with optional inter/intra/all filtering and residue-level annotation.

    Parameters:
        - chain_atoms_dict: {chain_id: [np.array([x,y,z]), ...]}
        - atom_residue_map: {chain_id: [(residue_id_str, resname, seqpos), ...]}
        - cutoff: distance threshold
        - mode: 'inter', 'intra', or 'all'
        - sequence_distance: minimum residue index separation for intra-chain
    """
    from scipy.spatial import KDTree

    contact_map = {}
    chain_ids = list(chain_atoms_dict.keys())
    total_pairs = (
        len(chain_ids) * (len(chain_ids) - 1) // 2 if mode == "inter"
        else sum(1 for _ in chain_ids) if mode == "intra"
        else len(chain_ids) * (len(chain_ids) + 1) // 2
    )
    progress_counter = 0

    for i in range(len(chain_ids)):
        for j in range(i, len(chain_ids)):
            c1, c2 = chain_ids[i], chain_ids[j]
            is_same_chain = c1 == c2

            # Filter mode
            if mode == "inter" and is_same_chain:
                continue
            if mode == "intra" and not is_same_chain:
                continue

            coords1 = np.array(chain_atoms_dict[c1])
            coords2 = np.array(chain_atoms_dict[c2])
            res1 = atom_residue_map[c1]
            res2 = atom_residue_map[c2]

            if coords1.size == 0 or coords2.size == 0:
                continue

            tree = KDTree(coords2)
            residue_contacts = set()

            for idx1, atom1 in enumerate(coords1):
                close_indices = tree.query_ball_point(atom1, r=cutoff)
                for idx2 in close_indices:
                    res1_id, _, seq1 = res1[idx1]
                    res2_id, _, seq2 = res2[idx2]

                    # For intra-chain, skip close-in-sequence contacts
                    if is_same_chain and abs(seq1 - seq2) < sequence_distance:
                        continue

                    residue_contacts.add((res1_id, res2_id))

            if residue_contacts:
                contact_map[(c1, c2)] = {
                    "count": len(residue_contacts),
                    "residue_pairs": sorted(residue_contacts)
                }

            progress_counter += 1
            if progress is not None:
                progress.value = 40 + int((progress_counter / total_pairs) * 40)
            if status_label is not None:
                status_label.value = f"Computing contacts... ({progress_counter}/{total_pairs})"

    return contact_map


def compute_contacts(chain_atoms, cutoff=5.0, progress=None, status_label=None):
    """Compute inter-chain contacts using KDTree for spatial efficiency."""
    contact_map = {}
    chain_ids = list(chain_atoms.keys())
    total_pairs = (len(chain_ids) * (len(chain_ids) - 1)) // 2
    progress_counter = 0

    for i in range(len(chain_ids)):
        for j in range(i + 1, len(chain_ids)):
            chain1, chain2 = chain_ids[i], chain_ids[j]
            coords1 = np.array(chain_atoms[chain1])
            coords2 = np.array(chain_atoms[chain2])

            # Skip empty chains
            if coords1.size == 0 or coords2.size == 0:
                continue
            if coords1.ndim != 2 or coords2.ndim != 2 or coords1.shape[1] != 3 or coords2.shape[1] != 3:
                continue

            # Build KDTree and query contacts
            tree = KDTree(coords2)
            num_contacts = sum(len(tree.query_ball_point(pt, r=cutoff)) > 0 for pt in coords1)

            if num_contacts > 0:
                contact_map[(chain1, chain2)] = num_contacts

            # Progress bar update
            progress_counter += 1
            if progress is not None:
                progress.value = 40 + int((progress_counter / total_pairs) * 40)
            if status_label is not None:
                status_label.value = f"Computing contacts... ({progress_counter}/{total_pairs})"

    return contact_map


def plot_and_save_heatmap(contact_map, output_file):
    """Plot and save a heatmap of contacts."""
    chains = sorted(set([c for pair in contact_map for c in pair]))
    matrix = np.zeros((len(chains), len(chains)))
    idx = {c: i for i, c in enumerate(chains)}

    for (c1, c2), val in contact_map.items():
      i, j = idx[c1], idx[c2]
      count = val["count"] if isinstance(val, dict) else val
      matrix[i, j] = matrix[j, i] = float(count)


    plt.figure(figsize=(8, 6))
    sns.heatmap(matrix, xticklabels=chains, yticklabels=chains, cmap="Blues", square=True)
    plt.xlabel("Chain ID")
    plt.ylabel("Chain ID")
    plt.title("Contact Heatmap")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Heatmap saved as {output_file}")


def get_chain_palette(n):
    """Return a distinct color palette with n colors."""
    import seaborn as sns
    return sns.color_palette("hls", n).as_hex()

def plot_and_save_interactive_network(contact_map, output_file, edge_style, cluster_map=None, weight_contacts=False, exclude_chains=None):
    """
    Create and save an interactive network visualization using ipysigma.

    Args:
        contact_map: dict with contact data
        output_file: str, path to save HTML
        edge_style: str, 'curve', 'line', etc.
        cluster_map: dict, optional, maps chain_id -> cluster_id (for coloring)
        weight_contacts: bool, if True, scale node size based on contact count
        exclude_chains: set of chain IDs to exclude from the graph
    """
    G = nx.Graph()
    exclude_chains = set(exclude_chains) if exclude_chains else set()

    chains = sorted({c for pair in contact_map for c in pair if c not in exclude_chains})
    color_keys = sorted({cluster_map.get(c, c) for c in chains}) if cluster_map else chains
    chain_palette = get_chain_palette(len(color_keys))
    palette = {k: chain_palette[i] for i, k in enumerate(color_keys)}

    for (chain1, chain2), contact_info in contact_map.items():
        if chain1 in exclude_chains or chain2 in exclude_chains:
            continue

        weight = contact_info["count"] if isinstance(contact_info, dict) else int(contact_info)
        node_size = weight if weight_contacts else 100

        color1 = cluster_map.get(chain1, chain1) if cluster_map else chain1
        color2 = cluster_map.get(chain2, chain2) if cluster_map else chain2

        G.add_node(chain1, color=color1, node_size=node_size)
        G.add_node(chain2, color=color2, node_size=node_size)

        if G.has_edge(chain1, chain2):
            G[chain1][chain2]['weight'] += weight
        else:
            G.add_edge(chain1, chain2, weight=weight)

    for _, data in G.nodes(data=True):
        data["node_size"] = int(data["node_size"])
    for _, _, data in G.edges(data=True):
        data["weight"] = int(data["weight"])

    Sigma.write_html(
        G,
        output_file,
        fullscreen=True,
        node_color="color",
        node_metrics=["louvain"],
        node_size="node_size",
        node_size_range=(10, 20),
        max_categorical_colors=30,
        edge_size="weight",
        edge_size_range=(5, 15),
        default_edge_type=edge_style,
        node_border_color_from="node",
        default_node_label_size=16,
        node_color_palette=palette
    )
    print(f"✅ Network saved: {output_file}")



def plot_networks_with_selfloops(contact_map, output_html, edge_style, output_tiff, dpi=300):
    """
    Generates two network visualizations:
    1. An interactive HTML file without self-loops using ipysigma.
    2. A high-resolution TIFF image with self-loops using NetworkX and Matplotlib.

    Parameters:
    - contact_map: dict
        Dictionary where keys are (chain1, chain2) and values are either a count or dict with "count".
    - output_html: str
        Path to the interactive HTML file.
    - output_tiff: str
        Path to the static .tiff image file.
    - dpi: int
        Resolution for the static image (default 300).
    """
    # --- Build graph ---
    G = nx.Graph()

    # Extract unique chain IDs
    chains = sorted(set([c for pair in contact_map for c in pair]))

    # Use a colorblind-friendly palette
    chain_palette = get_chain_palette(len(chains))
    palette = {chain: chain_palette[i] for i, chain in enumerate(chains)}

    # Add nodes and edges
    for (chain1, chain2), contact_info in contact_map.items():
        weight = contact_info["count"] if isinstance(contact_info, dict) else int(contact_info)
        G.add_node(chain1, color=chain1, node_size=weight)
        G.add_node(chain2, color=chain2, node_size=weight)

        if G.has_edge(chain1, chain2):
            G[chain1][chain2]['weight'] += weight
        else:
            G.add_edge(chain1, chain2, weight=weight)

    # Normalize attributes
    for _, data in G.nodes(data=True):
        data["node_size"] = int(data["node_size"])
    for _, _, data in G.edges(data=True):
        data["weight"] = int(data["weight"])

    # --- Interactive HTML without self-loops ---
    G_no_selfloops = G.copy()
    self_loops = list(nx.selfloop_edges(G_no_selfloops))
    G_no_selfloops.remove_edges_from(self_loops)

    Sigma.write_html(
        G,
        output_html,
        fullscreen=True,
        node_color="color",  # category name
        node_metrics=["louvain"],
        node_size="node_size",
        #node_size_range=(10, 20),
        max_categorical_colors=30,
        edge_size="weight",
        edge_size_range=(5, 15),
        default_edge_type=edge_style,
        node_border_color_from="node",
        default_node_label_size=16,
        node_color_palette=palette  # category -> color
    )
    print(f"✅ Interactive network saved as {output_html}")

    # --- Static .tiff with self-loops ---
    pos = nx.spring_layout(G, seed=42)  # Reproducible layout

    plt.figure(figsize=(10, 10), dpi=dpi)

    # Use same node color palette
    node_colors = [palette.get(node, "#999999") for node in G.nodes()]
    nx.draw_networkx_nodes(G, pos, node_size=300, node_color=node_colors)
    nx.draw_networkx_labels(G, pos, font_size=10)

    # Separate self-loops
    self_loop_edges = list(nx.selfloop_edges(G))
    other_edges = [edge for edge in G.edges() if edge not in self_loop_edges]

    # Draw inter-chain edges
    nx.draw_networkx_edges(G, pos, edgelist=other_edges, width=1.5, edge_color='gray')

    # Draw self-loops as red circles
    for node in G.nodes():
        if G.has_edge(node, node):
            loop_weight = G[node][node]['weight']
            loop = plt.Circle(pos[node], 0.05, color='red', fill=False, linewidth=1.5)
            plt.gca().add_patch(loop)
            # Optional: annotate loop count
            # plt.text(pos[node][0], pos[node][1] + 0.07, f"{loop_weight}", fontsize=8, ha='center')

    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_tiff, format='tiff', dpi=dpi)
    plt.close()
    print(f"✅ High-resolution network image saved as {output_tiff}")


def plot_residue_level_network(contact_map, output_file="network_residue_level.html", min_contacts=10):
    """
    Generates and saves a residue-level interactive contact network using ipysigma.

    Args:
        contact_map: dict from compute_contacts_with_residues
        output_file: HTML output path
        min_contacts: minimum number of edges for a residue to be styled normally
    """
    import networkx as nx
    import matplotlib.pyplot as plt
    from ipysigma import Sigma

    G = nx.Graph()

    # Step 1: Extract all chains
    all_chains = set()
    for (chain1, chain2), data in contact_map.items():
        for res1, res2 in data.get("residue_pairs", []):
            all_chains.add(res1.split(":")[0])
            all_chains.add(res2.split(":")[0])
    sorted_chains = sorted(all_chains)

    # Step 2: Build color palette
    chain_palette = get_chain_palette(len(sorted_chains))
    palette = {chain: chain_palette[i] for i, chain in enumerate(sorted_chains)}

    # Step 3: Build graph (nodes first, to track all)
    all_residues = set()
    for (_, _), data in contact_map.items():
        for res1, res2 in data.get("residue_pairs", []):
            all_residues.add(res1)
            all_residues.add(res2)

    for residue in all_residues:
        chain = residue.split(":")[0]
        G.add_node(residue, label=residue, color=chain, node_size=5)

    # Step 4: Add edges
    for (_, _), data in contact_map.items():
        for res1, res2 in data.get("residue_pairs", []):
            if G.has_edge(res1, res2):
                G[res1][res2]["weight"] += 1
            else:
                G.add_edge(res1, res2, weight=1)

    # Step 5: Update node styles based on degree
    for node in G.nodes():
        deg = G.degree(node)
        if deg < min_contacts:
            G.nodes[node]["color"] = "#cccccc"  # light gray
            G.nodes[node]["node_size"] = 2      # smaller
        else:
            G.nodes[node]["node_size"] = 5      # default size

    # Convert attributes to int
    for _, data in G.nodes(data=True):
        data["node_size"] = int(data["node_size"])
    for _, _, data in G.edges(data=True):
        data["weight"] = int(data["weight"])

    # Step 6: Save interactive graph
    Sigma.write_html(
        G,
        output_file,
        fullscreen=True,
        node_color="color",
        node_size="node_size",
        node_size_range=(2, 8),
        edge_size="weight",
        edge_size_range=(1, 5),
        default_edge_type="curve",
        node_border_color_from="node",
        default_node_label_size=12,
        node_color_palette=palette  # for valid color categories
    )

    print(f"✅ Residue-level network saved as {output_file}")



def analyze_structure(file_path, cutoff=8):
    """Pipeline: Parse file, compute contacts, generate visualizations."""
    if not os.path.exists(file_path):
        print(f"❌ File not found: {file_path}")
        return

    print(f"🔍 Analyzing {file_path} with cutoff {cutoff} Å ...")
    try:
        chain_atoms = parse_structure(file_path, atom_mode=atom_selector.value)
        contact_map = compute_contacts(chain_atoms, cutoff)

        if not contact_map:
            print("⚠️ No inter-chain contacts found.")
            return

        plot_and_save_interactive_network(contact_map, NETWORK_HTML_FILE)
        plot_and_save_heatmap(contact_map, HEATMAP_FILE)
    except Exception as e:
        print(f"❌ Error during analysis: {e}")


def export_contacts_table(contact_map, output_file, structure_path):
    """
    Export residue-residue contacts to TSV including residue names and sequence numbers.

    Args:
        contact_map: dict with 'residue_pairs' as (res1_id, res2_id)
        output_file: output path to save the TSV
        structure_path: path to the PDB or CIF file to extract residue names
    """


    # Load structure
    if structure_path.endswith(".pdb"):
        parser = PDBParser(QUIET=True)
    elif structure_path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format.")

    structure = parser.get_structure("structure", structure_path)

    # Create lookup: "A:1002" -> ("ALA", 1002)
    residue_lookup = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                chain_id = chain.id
                resseq = residue.id[1]
                resname = residue.get_resname()
                key = f"{chain_id}:{resseq}"
                residue_lookup[key] = (resname, resseq)

    # Write TSV
    with open(output_file, "w", newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["chain1", "res1", "resname1", "chain2", "res2", "resname2"])

        for (chain1, chain2), data in contact_map.items():
            if "residue_pairs" not in data:
                continue
            for res1, res2 in data["residue_pairs"]:
                resname1, resseq1 = residue_lookup.get(res1, ("UNK", -1))
                resname2, resseq2 = residue_lookup.get(res2, ("UNK", -1))
                chain_id1, _ = res1.split(":")
                chain_id2, _ = res2.split(":")
                writer.writerow([chain_id1, resseq1, resname1, chain_id2, resseq2, resname2])

    print(f"✅ Contact table saved to: {output_file}")


def download_pdb(pdb_id, output_path_base):
    """
    Downloads a PDB structure file in .pdb format. If not available, tries .cif format.

    Args:
        pdb_id (str): The 4-character PDB ID (e.g., "1TUP").
        output_path_base (str): Base path without extension.

    Returns:
        str or None: Path to the downloaded file, or None if both downloads failed.
    """
    import requests

    pdb_id = pdb_id.strip().upper()

    # Try .pdb format first
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(pdb_url)
        response.raise_for_status()
        output_path = f"{output_path_base}.pdb"
        with open(output_path, 'w') as f:
            f.write(response.text)
        print(f"✅ Successfully downloaded {pdb_id}.pdb")
        return output_path
    except requests.exceptions.RequestException:
        print(f"⚠️  Failed to download {pdb_id}.pdb. Trying .cif format...")

    # Try .cif format as fallback
    cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try:
        response = requests.get(cif_url)
        response.raise_for_status()
        output_path = f"{output_path_base}.cif"
        with open(output_path, 'w') as f:
            f.write(response.text)
        print(f"✅ Successfully downloaded {pdb_id}.cif")
        return output_path
    except requests.exceptions.RequestException:
        print(f"❌ Failed to download both .pdb and .cif for {pdb_id}")
        return None

def extract_chain_sequences(file_path, output_fasta):
    """
    Extract chain sequences from a PDB or CIF file and save them in FASTA format.

    Args:
        file_path (str): Path to structure file (.pdb or .cif)
        output_fasta (str): Output FASTA file
    Returns:
        Dict mapping chain ID -> sequence
    """
    if file_path.endswith(".pdb"):
        parser = PDB.PDBParser(QUIET=True)
    elif file_path.endswith(".cif"):
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format.")

    structure = parser.get_structure("structure", file_path)
    ppb = PPBuilder()
    chain_seqs = {}

    for model in structure:
        for chain in model:
            chain_id = chain.id
            peptides = ppb.build_peptides(chain)
            seq = ''.join(str(p.get_sequence()) for p in peptides)
            if seq:
                chain_seqs[chain_id] = seq

    # Save to FASTA
    records = [SeqRecord(Seq(seq), id=chain, description="") for chain, seq in chain_seqs.items()]
    SeqIO.write(records, output_fasta, "fasta")

    return chain_seqs

import os
import subprocess

def find_similar_chains(pdb_fasta, reference_fasta, tmp_dir="/tmp/mmseqs_filter", min_seq_id=0.6,force=False):
    
    if force and os.path.exists(tmp_dir):
        print(f"🧹 Removing previous MMseqs output in {tmp_dir}")
        shutil.rmtree(tmp_dir)
        
    os.makedirs(tmp_dir, exist_ok=True)

    query_db = os.path.join(tmp_dir, "queryDB")
    target_db = os.path.join(tmp_dir, "targetDB")
    result_db = os.path.join(tmp_dir, "resultDB")
    result_tsv = os.path.join(tmp_dir, "result.tsv")

    cmds = [
        f"mmseqs createdb {pdb_fasta} {query_db} -v 0",
        f"mmseqs createdb {reference_fasta} {target_db} -v 0",
        f"mmseqs search {query_db} {target_db} {result_db} {tmp_dir}/tmp --min-seq-id {min_seq_id} -c 0.8 -v 0",
        f"mmseqs convertalis {query_db} {target_db} {result_db} {result_tsv} -v 0"
    ]

    for cmd in cmds:
        print(f"💻 Running: {cmd}")
        subprocess.run(cmd, shell=True, check=True)

    # Parse hits
    similar_chains = set()
    with open(result_tsv) as f:
        for line in f:
            query_id, *_ = line.strip().split('\t')
            similar_chains.add(query_id)

    return similar_chains



def load_fasta_sequences(fasta_path):
    """Load sequences from a FASTA file."""
    return {str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}


def run_mmseqs_clustering(fasta_file, output_dir, min_identity=0.8, force=False):
    """
    Run MMseqs2 clustering on the input FASTA sequences.

    Args:
        fasta_file: Path to FASTA file
        output_dir: Path to MMseqs2 output directory
        min_identity: Minimum sequence identity for clustering
        force: If True, remove previous MMseqs2 output before running

    Returns:
        dict: chain_id -> cluster_id
    """
    import subprocess
    import os

    if force and os.path.exists(output_dir):
        print(f"🧹 Removing previous MMseqs output in {output_dir}")
        shutil.rmtree(output_dir)

    os.makedirs(output_dir, exist_ok=True)
    db = os.path.join(output_dir, "seqDB")
    tmp = os.path.join(output_dir, "tmp")
    result = os.path.join(output_dir, "clusterRes")

    cmds = [
        f"mmseqs createdb {fasta_file} {db} -v 0",
        f"mmseqs cluster {db} {result} {tmp} --min-seq-id {min_identity} -c 0.5 -v 0",
        f"mmseqs createtsv {db} {db} {result} {output_dir}/clusters.tsv -v 0"
    ]

    for cmd in cmds:
        print(f"💻 Running: {cmd}")
        subprocess.run(cmd, shell=True, check=True)

    # Parse clusters.tsv
    cluster_map = {}
    with open(f"{output_dir}/clusters.tsv") as f:
        for line in f:
            rep, member = line.strip().split('\t')
            cluster_map[member] = rep

    return cluster_map



import os
import time
import argparse

def main():
    parser = argparse.ArgumentParser(description="ProtCNet CLI (refactored from Colab)")

    # File input
    parser.add_argument("--pdb_id", type=str, help="PDB ID to download")
    parser.add_argument("--input_file", type=str, help="Path to local .pdb or .cif file")
    parser.add_argument("--exclude_fasta", type=str, help="FASTA file with sequences to exclude if similar (>0.8)")
    parser.add_argument( "--exclude_chain", nargs="+", default=[], help="List of chain IDs to exclude from contact network (e.g. A B C)")

    # Contact computation
    parser.add_argument("--cutoff", type=float, default=8.0, help="Distance cutoff (Å)")
    parser.add_argument("--mode", choices=["inter", "intra", "all"], default="inter", help="Contact type")
    parser.add_argument("--sequence_distance", type=int, default=15, help="Minimum residue separation for intra-chain contacts")

    # Atom selection and filtering
    parser.add_argument("--atom_mode", choices=["all", "ca", "cb"], default="ca", help="Atoms to use for contact calculation")
    parser.add_argument("--interaction_filter", choices=["all", "electrostatic", "hydrophobic"], default="all", help="Residue filtering")

    # Flags
    parser.add_argument("--track_residues", action="store_true", help="Track contacting residues")
    parser.add_argument("--cluster_chains", action="store_true", help="Cluster chains by sequence similarity (MMseqs2)")
    parser.add_argument("--weight_contacts", action="store_true", help="Scale node size by contact count")
    parser.add_argument("--residue_level_net", action="store_true", help="Plot residue-level contact network")
    parser.add_argument("--force", action="store_true", help="Force overun mmseqs2, will delete previous mmseqs2 run.")
    
    # Residue-level network
    parser.add_argument("--residue_level_cutoff", type=float, default=10.0, help="Minimum contacts for residue-level node display")
    
    # Output
    parser.add_argument("--network_file", type=str, default="network_visualization.html", help="Output file for network")
    parser.add_argument("--edge_style", choices=["rectangle", "line", "curve"], default="curve", help="Edge style in network")

    args = parser.parse_args()

    start_time = time.time()
    print("🔍 Starting analysis...")

    file_path = None
    if args.input_file:
        file_path = args.input_file
    elif args.pdb_id:
        file_path = f"/tmp/{args.pdb_id.upper()}"
        file_path = download_pdb(args.pdb_id, file_path)
        if file_path is None:
            print("❌ File download failed.")
            return
    else:
        print("❌ No input file provided.")
        return

    print("📦 Parsing structure...")

    exclude_chains = set()
    if args.exclude_fasta:
        print("🧬 Filtering chains matching sequences in FASTA...")
    
        pdb_fasta = "/tmp/pdb_chains.fasta"
        extract_chain_sequences(file_path, pdb_fasta)
    
        exclude_chains = find_similar_chains(
            pdb_fasta,
            args.exclude_fasta,
            tmp_dir="/tmp/mmseqs_filter",
            min_seq_id=0.6,
            force = args.force
        )

    if args.exclude_chain:
        exclude_chains.update(args.exclude_chain)
    
    for chain_id in exclude_chains:
        print(f"❌ Excluding chain {chain_id}.")

    try:
        if args.track_residues:
            chain_atoms, residue_map = parse_structure(
                file_path,
                atom_mode=args.atom_mode,
                residue_filter=args.interaction_filter,
                with_residues=True
            )
        else:
            chain_atoms = parse_structure(
                file_path,
                atom_mode=args.atom_mode,
                residue_filter=args.interaction_filter,
                with_residues=False
            )
    except Exception as e:
        print(f"❌ Parsing failed: {e}")
        return

    cluster_map = None
    if args.cluster_chains:
        print("🧬 Extracting chain sequences and clustering...")
        try:
            extract_chain_sequences(file_path, "/tmp/chains.fasta")
            cluster_map = run_mmseqs_clustering("/tmp/chains.fasta", "/tmp/mmseqs_clusters", 0.8, args.force)
        except Exception as e:
            print(f"❌ MMseqs2 clustering failed: {e}")
            return

    print("📊 Computing contacts...")
    
    if exclude_chains:
        chain_atoms = {k: v for k, v in chain_atoms.items() if k not in exclude_chains}
        if args.track_residues:
            residue_map = {k: v for k, v in residue_map.items() if k not in exclude_chains}

    try:
        if args.track_residues:
            contact_map = compute_contacts_with_residues(
                chain_atoms,
                residue_map,
                cutoff=args.cutoff,
                mode=args.mode,
                sequence_distance=args.sequence_distance
            )
        else:
            contact_map = compute_contacts(
                chain_atoms,
                cutoff=args.cutoff
            )
    except Exception as e:
        print(f"❌ Contact computation failed: {e}")
        return

    if not contact_map:
        print("⚠️ No contacts found.")
        return

    print("🌐 Generating interactive network...")
    try:
        if args.mode == 'inter':
            plot_and_save_interactive_network(
            contact_map,
            output_file=args.network_file,
            edge_style=args.edge_style,
            cluster_map=cluster_map,
            weight_contacts=args.weight_contacts,
            exclude_chains=exclude_chains
        )

        else:
            plot_networks_with_selfloops(
                contact_map,
                args.network_file,
                args.edge_style,
                args.network_file.replace('.html', '.tiff')
            )
    except Exception as e:
        print(f"❌ Network plot failed: {e}")
        return

    print("🧯 Generating heatmap...")
    try:
        plot_and_save_heatmap(contact_map, 'heatmap.tiff')
    except Exception as e:
        print(f"❌ Heatmap plot failed: {e}")
        return

    if args.track_residues:
        try:
            print("📑 Exporting contact table...")
            export_contacts_table(
                contact_map,
                args.network_file.replace(".html", "_contacts.tsv"),
                file_path
            )
            if args.residue_level_net:
                print("🧠 Generating residue-level network...")
                plot_residue_level_network(
                    contact_map,
                    output_file="network_residue_level.html",
                    min_contacts=args.residue_level_cutoff
                )
        except Exception as e:
            print(f"⚠️ Failed to export contact table or residue network: {e}")

    elapsed = time.time() - start_time
    print(f"✅ Analysis complete! Runtime: {elapsed:.2f} seconds")

if __name__ == "__main__":
    main()

