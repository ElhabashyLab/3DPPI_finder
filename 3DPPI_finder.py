# =========================================
# 3DPPI_finder: Structural PPI Detection
# =========================================
# Author: Hadeer Elhabashy
# Date: 2026
#
# Description:
# - Identifies suitable PDB templates for protein pairs using BLAST.
# - Detects experimentally resolved protein complexes.
# - Maps protein–protein interfaces based on Cα distance criteria.
# - Generates per-pair interface reports and a global summary.
# =========================================


#!/usr/bin/env python3

import os
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from Bio.PDB import MMCIFParser
from scipy.spatial.distance import cdist
import requests
import sys  # for simple progress bar

# ============================================================
# ====================== CONFIGURATION =======================
# ============================================================

# Input
PAIR_FILE = "/media/elhabashy/Elements/scripts/testPPI.csv"   # columns: uid1, uid2

# Base paths
BASE_DIR = Path("/media/elhabashy/Elements/scripts/test")
PROTEIN_DIR = BASE_DIR / "proteins"
RESULT_DIR = BASE_DIR / "3DPPI"

# BLAST
BLAST_BIN = Path.home() / "Documents/software/blast/ncbi-blast-2.16.0+/bin/blastp"
BLAST_DB = Path.home() / "Documents/software/blast/blastdb/pdbaa_01_10_2025/pdbaa"
EVALUE = "1e-5"
IDENTITY_CUTOFF = 50

# UniProt
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{}.fasta"

# RCSB
RCSB_CIF_URL = "https://files.rcsb.org/view/{}.cif"

# Interface detection
DISTANCE_CUTOFF = 10.0
MIN_CONTACTS = 10

# ============================================================
# ======================= HELPERS ============================
# ============================================================

def run_cmd(cmd):
    """Run a shell command with check."""
    subprocess.run(cmd, shell=True, check=True)

def download_fasta(uid, out_file):
    """Download UniProt FASTA file if not already present."""
    if out_file.exists():
        return
    #print(f"  ↳ Downloading FASTA for {uid}...")
    run_cmd(f"wget -q -O {out_file} {UNIPROT_FASTA_URL.format(uid)}")

def run_blast(uid, fasta, out_file):
    """Run BLAST against the local PDB sequence database."""
    if out_file.exists():
        return
    #print(f"  ↳ Running BLAST for {uid}...")
    cmd = (
        f"{BLAST_BIN} "
        f"-query {fasta} "
        f"-db {BLAST_DB} "
        f"-evalue {EVALUE} "
        f"-outfmt 6 "
        f"-out {out_file}"
    )
    run_cmd(cmd)

def parse_blast(blast_file):
    """Parse BLAST tabular output and filter by identity cutoff."""
    cols = [
        "query", "subject", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    df = pd.read_csv(blast_file, sep="\t", names=cols)
    df["pdb"] = df["subject"].str.split("_").str[0]
    df["chain"] = df["subject"].str.split("_").str[1]
    df = df.drop_duplicates("subject")
    df = df[df["pident"] >= IDENTITY_CUTOFF]
    return df

def download_cif(pdb_id, out_file):
    """Download PDB CIF file if not already present."""
    if out_file.exists():
        return
    print(f"  ↳ Downloading CIF for {pdb_id}...")
    url = RCSB_CIF_URL.format(pdb_id.upper())
    r = requests.get(url)
    r.raise_for_status()
    out_file.write_bytes(r.content)

def get_ca_coords(structure, chain_id):
    """Return CA coordinates of a given chain."""
    coords = []
    for model in structure:
        if chain_id not in model:
            continue
        chain = model[chain_id]
        for res in chain:
            if "CA" in res:
                coords.append(res["CA"].coord)
    return np.array(coords)

def detect_interface(coords1, coords2):
    """Detect interface between two sets of coordinates."""
    if len(coords1) == 0 or len(coords2) == 0:
        return False, 0
    dists = cdist(coords1, coords2)
    n_contacts = np.sum(dists <= DISTANCE_CUTOFF)
    return n_contacts >= MIN_CONTACTS, n_contacts

def print_progress(current, total, prefix=""):
    """Simple progress bar for loops."""
    bar_len = 40
    filled_len = int(round(bar_len * current / float(total)))
    percents = round(100.0 * current / float(total), 1)
    bar = "#" * filled_len + "-" * (bar_len - filled_len)
    sys.stdout.write(f"\r{prefix} [{bar}] {percents}% ({current}/{total})")
    sys.stdout.flush()
    if current == total:
        print()  # newline at end

# ============================================================
# ========================= MAIN =============================
# ============================================================

def main():

    print("🚀 Starting 3D PPI interface detection script...")
    PROTEIN_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    pairs = pd.read_csv(PAIR_FILE)
    uids = sorted(set(pairs["uid1"]) | set(pairs["uid2"]))

    print(f"🔹 Found {len(uids)} unique proteins, starting BLAST step...")

    blast_cache = {}

    # ---------------- BLAST ALL PROTEINS ----------------
    total_proteins = len(uids)
    for i, uid in enumerate(uids, 1):
        # print progress bar with UID
        bar_len = 40
        filled_len = int(round(bar_len * i / float(total_proteins)))
        percents = round(100.0 * i / float(total_proteins), 1)
        bar = "#" * filled_len + "-" * (bar_len - filled_len)
        sys.stdout.write(f"\r[BLAST] {uid} [{bar}] {percents}% ({i}/{total_proteins})")
        sys.stdout.flush()

        pdir = PROTEIN_DIR / uid
        pdir.mkdir(exist_ok=True)

        fasta = pdir / f"{uid}.fasta"
        blast_out = pdir / f"{uid}.blast.tsv"

        # print messages inline with progress bar
        sys.stdout.write(f"  Downloading FASTA...")
        sys.stdout.flush()
        download_fasta(uid, fasta)
        
    
        sys.stdout.write(f"  Running BLAST...")
        sys.stdout.flush()
        run_blast(uid, fasta, blast_out)
       

        blast_cache[uid] = parse_blast(blast_out)

    print("\n🔹 BLAST step completed. Starting structural PPI detection...")

    summary = []
    total_pairs = len(pairs)

    # ---------------- STRUCTURAL PPI + INTERFACE ----------------
    for idx, row in pairs.iterrows():
        uid1, uid2 = row["uid1"], row["uid2"]
        print_progress(idx + 1, total_pairs, prefix=f"[PPI] {uid1}–{uid2}")

        df1 = blast_cache[uid1]
        df2 = blast_cache[uid2]

        merged = df1.merge(df2, on="pdb", suffixes=("_A", "_B"))
        merged = merged[merged["subject_A"] != merged["subject_B"]]

        # remove symmetric duplicates
        merged["pair"] = merged.apply(
            lambda x: "_".join(sorted([str(x["subject_A"]), str(x["subject_B"])])),
            axis=1
        )
        merged = merged.drop_duplicates("pair")
        if merged.empty:
            continue

        # ---- Prepare directories ----
        complex_dir = RESULT_DIR / f"{uid1}_{uid2}"
        complex_pdb_dir = complex_dir / "complex_pdb"
        complex_pdb_dir.mkdir(parents=True, exist_ok=True)

        # ---- Process first PDB only ----
        pdb_id = merged.iloc[0]["pdb"]
        merged = merged[merged["pdb"] == pdb_id].copy()
        merged["no_pairs"] = len(merged)

        cif_file = complex_pdb_dir / f"{pdb_id}.cif"
        download_cif(pdb_id, cif_file)

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(pdb_id, cif_file)

        interfaces = []
        contacts = []

        for _, r in merged.iterrows():
            coordsA = get_ca_coords(structure, r["chain_A"])
            coordsB = get_ca_coords(structure, r["chain_B"])
            is_iface, n_cont = detect_interface(coordsA, coordsB)
            interfaces.append(is_iface)
            contacts.append(n_cont)

        merged["interface"] = interfaces
        merged["n_contacts"] = contacts

        # ---- Interface architecture ----
        iface = merged[merged["interface"]]
        if not iface.empty:
            arch_A = iface["chain_A"].nunique()
            arch_B = iface["chain_B"].nunique()
            interface_arch = f"{arch_A}:{arch_B}"
            chain1 = iface.iloc[0]["chain_A"]
            chain2 = iface.iloc[0]["chain_B"]
        else:
            interface_arch = "0:0"
            chain1 = ""
            chain2 = ""

        # ---- Save per-pair CSV ----
        out_csv = complex_dir / f"{uid1}-{uid2}_with_interface.csv"
        merged["chain1"] = merged["chain_A"]
        merged["chain2"] = merged["chain_B"]
        merged.to_csv(out_csv, index=False)

        # ---- Append to summary ----
        summary.append({
            "uid1": uid1,
            "uid2": uid2,
            "complex_pdb": pdb_id,
            "no_pairs": merged.iloc[0]["no_pairs"],
            "interface_pairs": iface.shape[0],
            "interface_architecture": interface_arch,
            "chain1": chain1,
            "chain2": chain2
        }) 

        # ---- Cleanup CIF and directory ----
        cif_file.unlink(missing_ok=True)
        if complex_pdb_dir.exists():
            for f in complex_pdb_dir.glob("*"):
                f.unlink(missing_ok=True)
            complex_pdb_dir.rmdir()

    # ---------------- SUMMARY ----------------
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(RESULT_DIR / "3Dppi_summary.csv", index=False)

    print("\n✅ 3D PPI interface detection complete!")
    print(f"Summary saved to: {RESULT_DIR / '3Dppi_summary.csv'}")

if __name__ == "__main__":
    main()
