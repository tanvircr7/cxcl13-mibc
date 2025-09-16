# CXCL13 as a Prognostic Biomarker in Muscle-Invasive Bladder Cancer (TCGA-BLCA)

Reproducible analysis and figures for our study evaluating CXCL13 expression as a prognostic factor in muscle-invasive bladder cancer (MIBC) using the TCGA-BLCA cohort. The repository provides a single-command R script that loads data, performs survival analyses (OS/PFS), fits multivariable Cox models, and exports publication-quality plots and sample lists.

## 🔬 Study Snapshot

| Parameter | Details |
|-----------|---------|
| **Hypothesis** | Higher CXCL13 expression associates with improved survival |
| **Cohort** | TCGA-BLCA |
| **Primary Analyses** | • 75th percentile split of CXCL13 expression (high vs. low)<br>• Kaplan–Meier (OS, PFS) with 5-year administrative censoring<br>• Cox proportional hazards (univariable + multivariable adjusted for age, sex, stage; Stage I/II vs III/IV collapse provided) |
| **Outputs** | KM curves with HR/CIs and risk tables, forest plots, and high/low sample lists for downstream GSEA |

## 📁 Repository Structure

```
.
├── CXCL13_expression.R                                    # main analysis script (one command)
├── data/
│   └── tcga_pancanceraltas_with_CXCL13subgroups.csv      # input CSV (or add instructions)
├── output/                                                # results written here (created automatically)
│   ├── km_os_CXCL13_annot.png
│   ├── km_dfs_CXCL13_annot.png
│   ├── forestplot_os_multivariable.png
│   ├── forestplot_pfs_multivariable.png
│   ├── CXCL13_highrisk_samples.csv
│   └── CXCL13_lowrisk_samples.csv
└── README.md
```

## 🧰 Requirements

- **R** ≥ 4.1 (tested with 4.3–4.5)
- **R packages**: `survival`, `survminer`, `ggplot2`, `ggpubr`
- **Linux users** may need system libraries for plotting & compilation

### Optional (Linux) – Faster Binary Installs

Add this at the top of your R session/script to use prebuilt binaries:

```r
# Ubuntu 22.04 (jammy)
options(repos = c(RSPM = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))

# Ubuntu 24.04 (noble)
# options(repos = c(RSPM = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))
```

**Install packages:**

```r
install.packages(c("survival","survminer","ggplot2","ggpubr"))
```

## 🚀 Quick Start

### 1. Clone the repository and enter the folder:

```bash
git clone <your-repo-url>.git
cd <your-repo-name>
```

### 2. Place data in `data/`:

```
data/tcga_pancanceraltas_with_CXCL13subgroups.csv
```

### 3. Run the analysis (creates `output/` automatically):

```bash
Rscript CXCL13_expression.R
```

### 4. Find results in `output/`:

- **KM plots**: `km_os_CXCL13_annot.png`, `km_dfs_CXCL13_annot.png`
- **Forest plots (multivariable)**: `forestplot_*`
- **Sample lists**: `CXCL13_highrisk_samples.csv`, `CXCL13_lowrisk_samples.csv`

## 📑 How to Cite

If you use this code or figures, please cite:


You may also include a license badge and DOI (Zenodo) once minted.

## 📄 License

Choose one and add `LICENSE` at repo root:

- **MIT** (code permissive)
- **CC BY 4.0** (data/figures attribution)

## ✉️ Contact


---

*Last updated: 2025*