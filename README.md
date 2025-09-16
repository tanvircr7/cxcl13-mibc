# CXCL13 as a Prognostic Biomarker in Muscle-Invasive Bladder Cancer (TCGA-BLCA)

Reproducible analysis and figures for our study evaluating **CXCL13 expression** as a prognostic factor in muscle-invasive bladder cancer (MIBC) using the TCGA-BLCA cohort. The repository provides a single-command R script that loads data, performs survival analyses (OS/PFS), fits multivariable Cox models, and exports publication-quality plots and sample lists.

---

## 🔬 Study Snapshot

| Parameter            | Details                                                                                                                                                                                                                                                       |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Hypothesis**       | Higher CXCL13 expression associates with improved survival                                                                                                                                                                                                    |
| **Cohort**           | TCGA-BLCA                                                                                                                                                                                                                                                     |
| **Primary Analyses** | • 75th percentile split of CXCL13 expression (high vs. low)<br>• Kaplan–Meier (OS, PFS) with 5-year administrative censoring<br>• Cox proportional hazards (univariable + multivariable adjusted for age, sex, stage; Stage I/II vs III/IV collapse provided) |
| **Outputs**          | KM curves with HR/CIs and risk tables, forest plots, and high/low sample lists for downstream GSEA                                                                                                                                                            |

---

## 📁 Repository Structure

```
.
├── CXCL13_expression.R                # main analysis script
├── data/
│   └── tcga_pancanceraltas_with_CXCL13subgroups.csv
├── CXCL13_highrisk_samples.csv        # high-expression cases for GSEA
├── CXCL13_lowrisk_samples.csv         # low-expression cases for GSEA
├── km_os_CXCL13_annot.png             # KM curve (OS)
├── km_dfs_CXCL13_annot.png            # KM curve (PFS)
├── km_os_CXCL13_risktable.png         # OS risk table
├── km_dfs_CXCL13_risktable.png        # PFS risk table
├── km_os_CXCL13_combined.png          # OS combined KM
├── km_dfs_CXCL13_combined.png         # PFS combined KM
├── Rplots.pdf                         # default R autosave (may be empty)
├── plots/
│   ├── forestplot_os_multivariable.png
│   ├── forestplot_pfs_multivariable.png
│   ├── km_os_CXCL13_annot.png         # duplicate copy under /plots
│   └── km_dfs_CXCL13_annot.png
└── README.md
```

---

## 🧰 Requirements

* **R** ≥ 4.1 (tested with 4.3–4.5)
* **R packages**: `survival`, `survminer`, `ggplot2`, `ggpubr`
* **Linux users** may need system libraries for plotting & compilation

### Optional (Linux) – Faster Binary Installs

```r
# Ubuntu 22.04 (jammy)
options(repos = c(RSPM = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))

# Ubuntu 24.04 (noble)
# options(repos = c(RSPM = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))
```

Install packages:

```r
install.packages(c("survival","survminer","ggplot2","ggpubr"))
```

---

## 🚀 Quick Start

### 1. Clone the repository and enter the folder:

```bash
git clone <your-repo-url>.git
cd cxcl13-mibc
```

### 2. Place data in `data/`:

```
data/tcga_pancanceraltas_with_CXCL13subgroups.csv
```

### 3. Run the analysis (creates outputs automatically):

```bash
Rscript CXCL13_expression.R
```

### 4. Results

* **KM plots (root)**:
  `km_os_CXCL13_annot.png`, `km_dfs_CXCL13_annot.png`, plus risk tables and combined plots
* **Forest plots (multivariable)**:
  in `plots/`
* **Sample lists**:
  `CXCL13_highrisk_samples.csv`, `CXCL13_lowrisk_samples.csv`

---

## 📑 How to Cite

If you use this code or figures, please cite this repository.

---

## 📄 License

Choose one and add `LICENSE` at repo root:

* **MIT** (for code)
* **CC BY 4.0** (for data/figures attribution)

---

## ✉️ Contact

For questions, reach out via GitHub issues or email.

---

*Last updated: 2025*

---

👉 Do you want me to make the **README “leaner”** (less verbose, geared to collaborators) or keep this more **publication-style detailed**?
