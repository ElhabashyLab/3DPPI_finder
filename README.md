# 3DPPI_finder
3DPPI_finder is a Python tool designed to analyze protein-protein interactions at the structural level. 
Given a list of protein pairs, the script first checks whether experimentally resolved 3D structures are available in the Protein Data Bank (PDB). 
For pairs with available structures, it identifies if they form protein complexes and maps the interfaces between interacting proteins.

## Requirements

- **Python ≥ 3.8**  
- **Python packages:** `pandas`, `numpy`, `biopython`, `scipy`, `requests`  
- **BLAST+** (tested with 2.16.0+)  
- **Local PDB sequence database** for BLAST (e.g., `pdbaa`)
- **Internet connection** (UniProt and RCSB PDB APIs are queried at runtime).


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


## **Configuration** 

Edit the top section of the script to configure paths and parameters:

# =======================
# Input / Output settings
# =======================

Path to the input CSV file containing protein pairs. The file must contain two columns: uid1 and uid2 (UniProt IDs).
PAIR_FILE = "/path/to/pairs.csv"

Base directory where all outputs and intermediate files will be stored.
BASE_DIR = Path("/path/to/output")

 Directory used to cache downloaded FASTA sequences and BLAST results for individual UniProt IDs.
PROTEIN_DIR = BASE_DIR / "proteins"

Directory where all per-pair and summary results will be written.
RESULT_DIR = BASE_DIR / "3DPPI"


# BLAST configuration


Path to the BLASTP executable (BLAST+).
BLAST_BIN = Path("/path/to/blastp")

Path to the local PDB sequence database used for BLAST searches (e.g., pdbaa or an equivalent PDB-derived FASTA database).
BLAST_DB = Path("/path/to/pdbaa")

E-value threshold used to filter BLAST hits. Lower values enforce stricter sequence similarity.
EVALUE = "1e-5"

- Minimum sequence identity (%) required for a BLAST hit to be considered as a structural template candidate.
-   IDENTITY_CUTOFF = 50

- Interface detection parameters
-   Maximum Cα–Cα distance (in Ångström) used to define a contact between two protein chains.
-   DISTANCE_CUTOFF = 10.0

- Minimum number of Cα–Cα contacts required to classify two proteins as interacting in a 3D complex.
-   MIN_CONTACTS = 10


## **Input Format**

The input CSV (PAIR_FILE) must contain two columns:  
uid1,uid2  
P12345,Q67890  
A11111,B22222  
C33333,D44444

uid1 and uid2 are UniProt IDs of proteins to test for interaction.


## **Output**

**1. Per-pair CSV**

Contains BLAST hits and interface detection for each protein pair:

<RESULT_DIR>/<uid1>_<uid2>/<uid1>-<uid2>_with_interface.csv 

Columns include:

- subject_A, subject_B: PDB chain IDs
- interface: True if interface detected
- n_contacts: Number of Cα contacts
- chain1, chain2: Chains forming the interface


**2. Summary CSV** 

Aggregates results across all pairs:

<RESULT_DIR>/3Dppi_summary.csv


Columns include:
- uid1, uid2
- complex_pdb: PDB ID used
- no_pairs: Number of candidate chain pairs
- interface_pairs: Number of chain pairs with interface
- interface_architecture: e.g., "2:1" (chains in uid1 : chains in uid2)
- chain1, chain2: Representative interface chains


## **Usage**

Run the script directly:
```bash
python3 3Dppi_finder.py
```

## **Notes**

- BLAST results and downloaded FASTA files are cached in <PROTEIN_DIR>/<uid>/.

- Only the first PDB structure per protein pair is considered for interface detection.
  
- Interface Detection Strategy:
  - A fast ("coarse"/"dirty") interface detection is used to maintain scalability.
  - Interfaces are defined based on Cα–Cα distances between protein chains.
  - Two proteins are considered interacting if:
  -   * At least MIN_CONTACTS (default = 10) Cα–Cα pairs
  -   * Are within DISTANCE_CUTOFF Å (default = 10.0 Å).


## **License**
This script is licensed under the MIT License. See the LICENSE file for more details.  

## Author
* Hadeer Elhabashy
