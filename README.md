# 3DPPI_finder
3DPPI_finder is a Python tool designed to analyze protein-protein interactions at the structural level. 
Given a list of protein pairs, the script first checks whether experimentally resolved 3D structures are available in the Protein Data Bank (PDB). 
For pairs with available structures, it identifies if they form protein complexes and maps the interfaces between interacting proteins.


---

## Requirements

- **Python ≥ 3.8**  
- **Python packages:** `pandas`, `numpy`, `biopython`, `scipy`, `requests`  
- **BLAST+** (tested with 2.16.0+)  
- **Local PDB sequence database** for BLAST (e.g., `pdbaa`)  

---

## Installation

### **Clone the repository**

```bash
git clone https://github.com/ElhabashyLab/3DPPI_finder.git
cd 3DPPI_finder
```

### **Install Python dependencies**
```bash
pip install pandas numpy biopython scipy requests
```

Ensure BLAST+ is installed and a PDB sequence database is available:
- https://www.ncbi.nlm.nih.gov/books/NBK569861/
- https://www.ncbi.nlm.nih.gov/books/NBK569850/

---

## **Configuration** 

Edit the top section of the script to configure paths and parameters:
- PAIR_FILE = "/path/to/pairs.csv"      # CSV with columns: uid1, uid2
- BASE_DIR = Path("/path/to/output")
- PROTEIN_DIR = BASE_DIR / "proteins"
- RESULT_DIR = BASE_DIR / "3DPPI"

- BLAST_BIN = Path("/path/to/blastp")
- BLAST_DB = Path("/path/to/pdbaa")
- EVALUE = "1e-5"
- IDENTITY_CUTOFF = 50

- DISTANCE_CUTOFF = 10.0
- MIN_CONTACTS = 10


## **Input Format**

The input CSV (PAIR_FILE) must contain two columns:

uid1,	uid2
P12345,	Q67890
...,	...

uid1 and uid2 are UniProt IDs of proteins to test for interaction.


## **Output**

**1. Per-pair CSV**

Contains BLAST hits and interface detection for each protein pair:

<RESULT_DIR>/<uid1>_<uid2>/<uid1>-<uid2>_with_interface.csv 


Columns include:

*subject_A*, subject_B: PDB chain IDs

interface: True if interface detected

n_contacts: Number of Cα contacts

chain1, chain2: Chains forming the interface


**2. Summary CSV** 

Aggregates results across all pairs:

<RESULT_DIR>/3Dppi_summary.csv


Columns include:

uid1, uid2

complex_pdb: PDB ID used

no_pairs: Number of candidate chain pairs

interface_pairs: Number of chain pairs with interface

interface_architecture: e.g., "2:1" (chains in uid1 : chains in uid2)

chain1, chain2: Representative interface chains


## **Usage**

Run the script directly:
```bash
python3 3Dppi_finder.py
```

## **Notes**

- BLAST results and downloaded FASTA files are cached in <PROTEIN_DIR>/<uid>/.

- Only the first PDB structure per protein pair is considered for interface detection.

- The interface is defined by Cα atoms within DISTANCE_CUTOFF Å, requiring at least MIN_CONTACTS contacts.


**License**

MIT License. See LICENSE
 for details.
