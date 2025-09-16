# TCGA Clinical Data & Gene Expression Pipeline

Pipeline to clean TCGA clinical data and add gene expression subgroups.

## Requirements
```bash
pip install pandas numpy matplotlib seaborn
```

## Input Files
- Expression matrix: `data_mrna_seq_v2_rsem.txt` (genes as rows, samples as columns)
- Clinical data: TCGA PanCancer Atlas CSV from cBioPortal

## Usage

### Step 1: Generate Gene Expression Subgroups
Edit `generic_gener.py`:
```python
INPUT_FILE = "data_mrna_seq_v2_rsem.txt"
TARGET_GENE = "CXCL13"  # or "CD274", "MRC1", etc.
OUTPUT_CSV = "blca_tcga_CXCL13_enrichment.csv"
```

Run:
```bash
python generic_gener.py
```

### Step 2: Merge Clinical Data + Subgroups
Edit `cleanerCLINSHEET_subgroup.py`:
```python
clinical = "tcga_pancanceratlas_clinical.csv"
subgroups = "blca_tcga_CXCL13_enrichment.csv"  # from Step 1
output_path = "tcga_pancanceratlas_with_CXCL13subgroups.csv"
```

Run:
```bash
python cleanerCLINSHEET_subgroup.py
```

## What It Does
- **Step 1**: Creates gene expression subgroups (enrinched/not_enriched)
- **Step 2**: Cleans clinical data and merges with subgroups

### Cleaning Operations:
- Cleans column names (whitespace → _, strips non-alphanumerics)
- Selects key clinical columns (Age, Stage, OS, DFS, MSI/TMB, etc.)
- Standardizes IDs by trimming Sample_ID to 15 chars and merging on sample_id
- Writes a merged CSV that contains clinical variables + your subgroup labels
- Optionally, renders a boxplot (default shows <GENE> expression by subgroup)

## Output
Final CSV with cleaned clinical variables + gene expression subgroup labels.