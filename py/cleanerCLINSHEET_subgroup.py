import pandas as pd
from pathlib import Path

# --- Config (relative, portable) ---
CLIN_PATIENT = Path("data/raw/data_clinical_patient.txt")
CLIN_SAMPLE  = Path("data/raw/data_clinical_sample.txt")
LABELS_CSV   = Path("data/processed/blca_tcga_CXCL13_enrichment.csv")
OUT_CSV      = Path("data/processed/tcga_pancanceratlas_with_CXCL13subgroups.csv")

OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

def clean_cols(cols):
    out = []
    for c in cols:
        c = c.strip().replace(" ", "_")
        c = c.encode("ascii", "ignore").decode()
        c = c.replace("-", "_")
        out.append(c)
    return out

# Load cBioPortal TSVs (skip comment lines)
pat = pd.read_csv(CLIN_PATIENT, sep="\t", dtype=str, comment="#")
smp = pd.read_csv(CLIN_SAMPLE,  sep="\t", dtype=str, comment="#")

# Clean column names
pat.columns = clean_cols(pat.columns)
smp.columns = clean_cols(smp.columns)

# Merge patient + sample (Patient_ID is the usual key)
clin = smp.merge(pat, on="Patient_ID", how="left", suffixes=("_sample", "_patient"))

# Keep common analysis columns if present
keep = [
    "Patient_ID", "Sample_ID", "Diagnosis_Age",
    "Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code",
    "Overall_Survival_Months", "Overall_Survival_Status",
    "Disease_Free_Months", "Disease_Free_Status",
    "Radiation_Therapy",
    "Neoadjuvant_Therapy_Type_Administered_Prior_To_Resection_Text",
    "TMB_nonsynonymous", "MSIsensor_Score", "MSI_MANTIS_Score",
    "Aneuploidy_Score", "Mutation_Count", "Sex",
    "Race_Category", "Ethnicity_Category",
    "Progress_Free_Survival_Months", "Progression_Free_Status",
    "Fraction_Genome_Altered"
]
keep = [k for k in keep if k in clin.columns]
clin = clin[keep].copy()

# Merge labels
labels = pd.read_csv(LABELS_CSV)
clin["Sample_ID"] = clin["Sample_ID"].astype(str).str.strip().str[:15]
labels["sample_id"] = labels["sample_id"].astype(str).str.strip().str[:15]

merged = clin.merge(labels, left_on="Sample_ID", right_on="sample_id", how="left") \
             .drop(columns=["sample_id"])

merged.to_csv(OUT_CSV, index=False)
print(f"✅ Saved merged clinical + labels -> {OUT_CSV.resolve()}")
print(f"   Rows: {len(merged)}  with label cols: {[c for c in merged.columns if c.endswith('_status')]}")
