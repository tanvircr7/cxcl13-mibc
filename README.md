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

> **Note**: If you cannot redistribute the full CSV, keep a small example snippet in `data/` and add download instructions (see [Data](#data) section below).

## 🧰 Requirements

- **R** ≥ 4.1 (tested with 4.3–4.5)
- **R packages**: `survival`, `survminer`, `ggplot2`, `ggpubr`
- **Linux users** may need system libraries for plotting & compilation (see [Troubleshooting](#troubleshooting))

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

## 📊 What the Script Does

The analysis pipeline performs the following steps:

1. **Data Loading & Preprocessing**
   - Loads the CSV and constructs survival objects with 5-year (60-month) administrative censoring
   - Splits CXCL13 expression at the 75th percentile to define high vs low groups

2. **Survival Analysis**
   - Produces Kaplan–Meier curves with log-rank p-values, confidence bands, and HR/C-index annotations
   - Fits Cox models for both OS and PFS:
     - **Univariable** (CXCL13 high vs low)
     - **Multivariable** adjusted for age, sex, and stage
     - Includes a Stage I/II vs III/IV collapse to reduce separation issues

3. **Output Generation**
   - Generates forest plots with fallbacks if factor separation causes instability
   - Exports GSEA-ready sample lists for high/low CXCL13 groups

## 🗂 Data

The analysis expects a CSV with at least these columns (TCGA-style):

| Column | Description |
|--------|-------------|
| `Overall_Survival_Months` | Overall survival time |
| `Overall_Survival_Status` | Overall survival status |
| `Progress_Free_Survival_Months` | Progression-free survival time |
| `Progression_Free_Status` | Progression-free survival status* |
| `Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code` | Cancer stage |
| `Diagnosis_Age` | Age at diagnosis |
| `Sex` | Patient sex |
| `Sample_ID` | Sample identifier |
| `CXCL13_expression` | CXCL13 expression values |

> ***Note**: May also be named `Disease_Free_*` in some datasets

**If you cannot include TCGA data in the repo**, add a short `data/README.md` describing how to obtain it from GDC/Firehose and how to preprocess to this schema.

### Error: "no package called 'survminer'"

**Solution**: Install required packages from R:

```r
install.packages(c("survival","survminer","ggplot2","ggpubr"))
```

### Linux compile errors (e.g., png/jpeg/xml2, Fortran, nloptr)

Install system dependencies, then reinstall R packages.

**Ubuntu/Debian:**

```bash
sudo apt update
sudo apt install -y \
  build-essential gfortran pkg-config \
  libcurl4-openssl-dev libxml2-dev \
  libpng-dev libjpeg-dev libtiff5-dev \
  libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
  libblas-dev liblapack-dev libnlopt-dev
```

**Fedora/RHEL:**

```bash
sudo dnf groupinstall -y "Development Tools"
sudo dnf install -y \
  gcc-gfortran pkgconfig \
  libcurl-devel libxml2-devel \
  libpng-devel libjpeg-turbo-devel libtiff-devel \
  freetype-devel fontconfig-devel harfbuzz-devel fribidi-devel \
  blas-devel lapack-devel nlopt-devel
```

### KM plots warn about "colour : 'CXCL13 group'"

This is a harmless message from ggpubr/survminer; plots still render correctly.

### Cox model warnings: "coefficient may be infinite"

This is typical with sparse stage levels. Use the Stage2 collapsed plots included in outputs; they mitigate separation issues.

## 📑 How to Cite

If you use this code or figures, please cite:

> Hasan MT, et al. CXCL13 expression is associated with improved survival in muscle-invasive bladder cancer (TCGA-BLCA). (Conference presentation, 2025). GitHub repository: `<add your repo URL>`.

You may also include a license badge and DOI (Zenodo) once minted.

## 📄 License

Choose one and add `LICENSE` at repo root:

- **MIT** (code permissive)
- **CC BY 4.0** (data/figures attribution)

## ✉️ Contact

For questions or issues, please open a GitHub issue or contact:

**Mohammad Tanvir Hasan**  
📧 [mthasan@ualr.edu](mailto:mthasan@ualr.edu)

---

*Last updated: 2025*