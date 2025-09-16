import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# --- Config (relative, portable) ---
INPUT_FILE = Path("data/raw/data_mrna_seq_v2_rsem.txt")
TARGET_GENE = "CXCL13"              # change if you want another gene
QUANTILE = 0.75                     # 75th percentile split
OUTPUT_CSV = Path("data/processed/blca_tcga_CXCL13_enrichment.csv")

# --- I/O prep ---
OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

# --- Load expression (genes x samples) ---
df = pd.read_csv(INPUT_FILE, sep="\t", index_col=0)

if TARGET_GENE not in df.index:
    raise SystemExit(f"❌ Gene {TARGET_GENE} not found in {INPUT_FILE}")

# numeric expression for TARGET_GENE
gene_expr = pd.to_numeric(df.loc[TARGET_GENE], errors="coerce")
# standardize sample IDs to 15 chars (TCGA style)
gene_expr.index = gene_expr.index.str.slice(0, 15)

# threshold & labels
thr = gene_expr.quantile(QUANTILE)
labels = pd.DataFrame({
    "sample_id": gene_expr.index,
    f"{TARGET_GENE}_expression": gene_expr.values,
    f"{TARGET_GENE}_status": np.where(gene_expr >= thr, "enriched", "not_enriched"),
})

# save
labels.to_csv(OUTPUT_CSV, index=False)
print(f"✅ Saved labels -> {OUTPUT_CSV.resolve()}")
print(f"   Gene={TARGET_GENE}, quantile={QUANTILE:.2f}, threshold={thr:.4g}")
