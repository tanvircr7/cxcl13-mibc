#!/usr/bin/env Rscript

# --- 1) Packages -----------------------------------------------------------
if (!requireNamespace("survival", quietly=TRUE)) install.packages("survival")
if (!requireNamespace("survminer", quietly=TRUE)) install.packages("survminer")
library(survival)
library(survminer)

# --- 2) Load & preprocess --------------------------------------------------
df <- read.csv("./data/processed/tcga_pancanceraltas_with_CXCL13subgroups.csv",
               stringsAsFactors = FALSE)

# --- Optional: Filter by tumor stage --------------------------------------
# df <- df[df$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code %in% c("STAGE III","STAGE IV"), ]

# --- 3) Survival variables + 5-year (60 mo) administrative censoring ------
df$OS_time   <- suppressWarnings(as.numeric(df$Overall_Survival_Months))
df$OS_event  <- as.numeric(grepl("^1", df$Overall_Survival_Status))
df$OS_time60  <- ifelse(!is.na(df$OS_time) & df$OS_time > 60, 60, df$OS_time)
df$OS_event60 <- ifelse(is.na(df$OS_time), NA, ifelse(df$OS_time > 60, 0, df$OS_event))

df$DFS_time   <- suppressWarnings(as.numeric(df$Progress_Free_Survival_Months))
df$DFS_event  <- as.numeric(grepl("^1", df$Progression_Free_Status))
df$DFS_time60  <- ifelse(!is.na(df$DFS_time) & df$DFS_time > 60, 60, df$DFS_time)
df$DFS_event60 <- ifelse(is.na(df$DFS_time), NA, ifelse(df$DFS_time > 60, 0, df$DFS_event))

# --- 4) CXCL13 expression & groups ----------------------------------------
df$CXCL13_expression <- suppressWarnings(as.numeric(df$CXCL13_expression))
df$CXCL13_zscore     <- scale(df$CXCL13_expression)
df$CXCL13_status     <- ifelse(df$CXCL13_zscore > 1, "overexpressed", "normal")

# 70th percentile split (change 0.7 -> 0.75 for 75th)
threshold_70 <- quantile(df$CXCL13_expression, 0.7, na.rm = TRUE)
df$CXCL13    <- ifelse(df$CXCL13_expression >= threshold_70, "high", "low")
df$CXCL13    <- factor(df$CXCL13, levels = c("low","high"))

# Provide a robust alias 'CXCL13_group' in case other code uses it
df$CXCL13_group <- df$CXCL13

# --- 5) KM plots + HR/CI/C-index annotations ------------------------------
## OS KM
fit_os <- survfit(Surv(OS_time60, OS_event60) ~ CXCL13, data = df)
xmax_os <- max(df$OS_time60, na.rm = TRUE)
p_os <- ggsurvplot(
  fit_os, data = df,
  pval = TRUE, pval.coord = c(0.65 * xmax_os, 0.32),
  conf.int = TRUE,
  title = "OS by CXCL13-exp (5-year)",
  xlab = "Time (months)", ylab = "OS Probability",
  legend.title = "CXCL13 group",
  risk.table = TRUE,
  xlim = c(0, xmax_os)
)

## DFS KM
fit_dfs <- survfit(Surv(DFS_time60, DFS_event60) ~ CXCL13, data = df)
xmax_dfs <- 60
p_dfs <- ggsurvplot(
  fit_dfs, data = df,
  pval = TRUE, pval.coord = c(0.65 * xmax_dfs, 0.32),
  conf.int = TRUE,
  xlim = c(0, xmax_dfs),
  title = "PFS by CXCL13-exp (5-year)",
  xlab = "Time (months)", ylab = "PFS Probability",
  legend.title = "CXCL13 group",
  risk.table = TRUE
)

# Cox models for annotation (binary split matching KM)
cox_os  <- coxph(Surv(OS_time60, OS_event60) ~ CXCL13, data = df)
cox_dfs <- coxph(Surv(DFS_time60, DFS_event60) ~ CXCL13, data = df)
s_os  <- summary(cox_os)
s_dfs <- summary(cox_dfs)

# Extract HR, CI, and C-index
hr_os   <- s_os$conf.int[1,"exp(coef)"]; lcl_os <- s_os$conf.int[1,"lower .95"]; ucl_os <- s_os$conf.int[1,"upper .95"]
cidx_os <- s_os$concordance[1]
hr_dfs  <- s_dfs$conf.int[1,"exp(coef)"]; lcl_dfs <- s_dfs$conf.int[1,"lower .95"]; ucl_dfs <- s_dfs$conf.int[1,"upper .95"]
cidx_dfs <- s_dfs$concordance[1]

# Add annotations (place just under the p-value)
lab_os  <- sprintf("HR = %.2f (%.2f–%.2f)\nC-index = %.3f", hr_os,  lcl_os,  ucl_os,  cidx_os)
lab_dfs <- sprintf("HR = %.2f (%.2f–%.2f)\nC-index = %.3f", hr_dfs, lcl_dfs, ucl_dfs, cidx_dfs)
p_os$plot  <- p_os$plot  + annotate("text", x = 0.65 * xmax_os,  y = 0.21, label = lab_os,  hjust = 0, size = 4.8)
p_dfs$plot <- p_dfs$plot + annotate("text", x = 0.65 * xmax_dfs, y = 0.21, label = lab_dfs, hjust = 0, size = 4.8)
png("plots/km_os_CXCL13_annot.png",  width = 1800, height = 1500, res = 200)
print(p_os)   # draws the combined grob to the device
dev.off()

png("plots/km_dfs_CXCL13_annot.png", width = 1800, height = 1500, res = 200)
print(p_dfs)
dev.off()


# --- 6) Save plots (save actual ggplot objects) ----------------------------
ggsave("km_os_CXCL13_annot.png",  plot = p_os$plot,  width = 7, height = 6, dpi = 150)
ggsave("km_dfs_CXCL13_annot.png", plot = p_dfs$plot, width = 7, height = 6, dpi = 150)
# Optional: risk tables & combined figures
ggsave("km_os_CXCL13_risktable.png",  plot = p_os$table,  width = 7, height = 2.5, dpi = 150)
ggsave("km_dfs_CXCL13_risktable.png", plot = p_dfs$table, width = 7, height = 2.5, dpi = 150)
os_combined  <- arrange_ggsurvplots(list(p_os),  print = FALSE)
dfs_combined <- arrange_ggsurvplots(list(p_dfs), print = FALSE)
ggsave("km_os_CXCL13_combined.png",  plot = os_combined,  width = 7, height = 8, dpi = 150)
ggsave("km_dfs_CXCL13_combined.png", plot = dfs_combined, width = 7, height = 8, dpi = 150)

# --- 7) Multivariable Cox (Age + Sex + Stage adjusted) --------------------
# Encode covariates
df$Sex   <- factor(df$Sex)
df$Stage <- factor(df$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code)

keep_os  <- complete.cases(df[, c("OS_time60","OS_event60","CXCL13_group","Diagnosis_Age","Sex","Stage")])
keep_dfs <- complete.cases(df[, c("DFS_time60","DFS_event60","CXCL13_group","Diagnosis_Age","Sex","Stage")])

cox_multi_os  <- coxph(Surv(OS_time60,  OS_event60)  ~ CXCL13_group + Diagnosis_Age + Sex + Stage, data = df[keep_os, ])
cox_multi_dfs <- coxph(Surv(DFS_time60, DFS_event60) ~ CXCL13_group + Diagnosis_Age + Sex + Stage, data = df[keep_dfs, ])

cat("\nCoxPH OS (multivariable):\n");  print(summary(cox_multi_os))
cat("\nCoxPH PFS (multivariable):\n"); print(summary(cox_multi_dfs))

# Tidy one-line tables
fmt <- function(fit) {
  s <- summary(fit)
  data.frame(
    term = rownames(s$coefficients),
    HR   = round(s$conf.int[,"exp(coef)"], 2),
    LCL  = round(s$conf.int[,"lower .95"], 2),
    UCL  = round(s$conf.int[,"upper .95"], 2),
    p    = signif(s$coefficients[,"Pr(>|z|)"], 3),
    Cidx = round(s$concordance[1], 3),
    row.names = NULL
  )
}
cat("\nTidy OS:\n");  print(fmt(cox_multi_os))
cat("\nTidy PFS:\n"); print(fmt(cox_multi_dfs))

# --- 8) Export patient lists for GSEA --------------------------------------
if (!("Sample_ID" %in% colnames(df))) stop("❌ Column 'Sample_ID' not found in the dataframe.")
write.csv(df[df$CXCL13 == "high", c("Sample_ID", "CXCL13_expression")],
          "CXCL13_highrisk_samples.csv",
          row.names = FALSE)
write.csv(df[df$CXCL13 == "low",  c("Sample_ID", "CXCL13_expression")],
          "CXCL13_lowrisk_samples.csv",
          row.names = FALSE)
cat("\n✅ Saved patient lists for GSEA:\n - CXCL13_highrisk_samples.csv\n - CXCL13_lowrisk_samples.csv\n")

# --- Forest plots ----------------------------------------------------------
# --- Stage collapse: I/II vs III/IV ---------------------------------------
stage_raw <- toupper(trimws(as.character(
  df$Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code
)))
df$Stage2 <- ifelse(grepl("^STAGE[[:space:]]*I\\b",   stage_raw) |
                      grepl("^STAGE[[:space:]]*II\\b",  stage_raw), "I/II",
                    ifelse(grepl("^STAGE[[:space:]]*III\\b", stage_raw) |
                             grepl("^STAGE[[:space:]]*IV\\b",  stage_raw), "III/IV", NA))
df$Stage2 <- factor(df$Stage2, levels = c("I/II","III/IV"))

# Standardize Sex a bit (optional but helps avoid rare/unknown levels)
df$Sex <- trimws(df$Sex)
df$Sex[df$Sex %in% c("M","Male","male")] <- "Male"
df$Sex[df$Sex %in% c("F","Female","female")] <- "Female"
df$Sex[!(df$Sex %in% c("Male","Female"))] <- NA
df$Sex <- factor(df$Sex)

# --- Rebuild complete-case masks using Stage2 ------------------------------
keep_os  <- complete.cases(df[, c("OS_time60","OS_event60","CXCL13_group","Diagnosis_Age","Sex","Stage2")])
keep_dfs <- complete.cases(df[, c("DFS_time60","DFS_event60","CXCL13_group","Diagnosis_Age","Sex","Stage2")])

d_os  <- droplevels(df[keep_os,  c("OS_time60","OS_event60","CXCL13_group","Diagnosis_Age","Sex","Stage2")])
d_dfs <- droplevels(df[keep_dfs, c("DFS_time60","DFS_event60","CXCL13_group","Diagnosis_Age","Sex","Stage2")])

# (Safety) Drop any factor levels that still have pure separation
drop_sep_levels <- function(dat, time, event, fac) {
  if (!is.factor(dat[[fac]])) return(dat)
  tb <- with(dat, table(dat[[fac]], dat[[event]]))
  bad <- rownames(tb)[rowSums(tb == 0) > 0]  # any level with a zero column (all events or all censored)
  if (length(bad)) {
    message(sprintf("⚠️ Dropping %s level(s) with separation: %s", fac, paste(bad, collapse=", ")))
    dat <- dat[!(dat[[fac]] %in% bad), , drop = FALSE]
    dat <- droplevels(dat)
  }
  dat
}
d_os  <- drop_sep_levels(d_os,  "OS_time60",  "OS_event60",  "Sex")
d_os  <- drop_sep_levels(d_os,  "OS_time60",  "OS_event60",  "Stage2")
d_dfs <- drop_sep_levels(d_dfs, "DFS_time60", "DFS_event60", "Sex")
d_dfs <- drop_sep_levels(d_dfs, "DFS_time60", "DFS_event60", "Stage2")

# --- Refit multivariable Cox with Stage2 -----------------------------------
cox_multi_os  <- coxph(Surv(OS_time60,  OS_event60)  ~ CXCL13_group + Diagnosis_Age + Sex + Stage2, data = d_os)
cox_multi_dfs <- coxph(Surv(DFS_time60, DFS_event60) ~ CXCL13_group + Diagnosis_Age + Sex + Stage2, data = d_dfs)

cat("\nCoxPH OS (multivariable, Stage2):\n");  print(summary(cox_multi_os))
cat("\nCoxPH PFS (multivariable, Stage2):\n"); print(summary(cox_multi_dfs))

# --- Forest plots ----------------------------------------------------------
fp_os <- ggforest(
  cox_multi_os, data = d_os,
  main = "Hazard ratios for OS (multivariable, Stage I/II vs III/IV)",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1.0,
  refLabel = "Reference",
  noDigits = 2
)
print(fp_os)
ggsave("plots/forestplot_os_multivariable.png", plot = fp_os, width = 7, height = 6, dpi = 150)

fp_dfs <- ggforest(
  cox_multi_dfs, data = d_dfs,
  main = "Hazard ratios for PFS (multivariable, Stage I/II vs III/IV)",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1.0,
  refLabel = "Reference",
  noDigits = 2
)
print(fp_dfs)
ggsave("plots/forestplot_pfs_multivariable.png", plot = fp_dfs, width = 7, height = 6, dpi = 150)

