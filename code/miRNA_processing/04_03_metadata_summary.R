rm(list = ls())

# Usage (run from project root):
# Rscript code/_miRNA_processing/04_03_metadata_summary.R 11
#
# Produces cohort characteristics (n, age, sex, APOE4, PMI, RIN, pathology) and
# summary statistics of the main findings (thresholds, effect sizes) needed to
# answer the reviewer request for the Abstract. Outputs go to Output<setup>/metadata/.

args <- commandArgs(trailingOnly = TRUE)
setup <- if (length(args) >= 1 && !is.na(args[1])) args[1] else ""

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(gtsummary)
})

data.dir   <- "data/preprocessed/"
demir.dir  <- "output/deMirs/"          # written by 02_deMirs.R and 03_01_miRNA_PathologyOddsModel.R
output.dir <- paste0("output", setup, "/")
meta.dir   <- file.path(output.dir, "metadata")
if (!dir.exists(meta.dir)) dir.create(meta.dir, recursive = TRUE, showWarnings = FALSE)

fdr.cut <- 0.05

covs <- as.data.frame(readRDS(file.path(data.dir, "mirClinical.RDS")))




#Defining resilience (identical to 02_deMirs.R)

cogdx.resils    <- c(1)
cogdx.AD        <- c(4)
cogdx.normal    <- c(1)

braak.resils    <- c(3,4,5,6)
braak.AD        <- c(3,4,5,6)

plq.resils      <- c(1,2)
plq.AD          <- c(1,2)
age.cut         <- 80
age.cut.AD      <- 80
age.cut.norm    <- 0


covs$age_death <- suppressWarnings(as.numeric(covs$age_death))

covs   <- mutate(covs, RESIL = "OTHER")

resilients     <- covs %>% dplyr::filter(.,
                                                 cogdx %in% cogdx.resils,
                                                 braaksc %in% braak.resils,
                                                 ceradsc %in% plq.resils,
                                                 age_death >= age.cut) %>%
  mutate(.,RESIL = "Resilient")#%>% nrow

AD     <- covs %>% filter(.,
                                          cogdx %in% cogdx.AD,
                                          braaksc %in% braak.AD,
                                          ceradsc %in% plq.AD,
                                          age_death >= age.cut.AD) %>%
  mutate(.,RESIL = "AD")

normals     <- covs %>% filter(.,
                                       cogdx %in% cogdx.normal,
                                       braaksc %in% c(0,1,2),
                                       ceradsc %in% c(3,4),
                                       age_death >= age.cut.norm) %>%
  mutate(.,RESIL = "Control")

presym     <- covs %>% filter(.,
                                      cogdx %in% c(2),
                                      braaksc %in% braak.resils,
                                      ceradsc %in% plq.resils) %>%
  mutate(.,RESIL = "Presym")

covs$RESIL[covs$mirna_id %in% resilients$mirna_id] <- "Resilient"
covs$RESIL[covs$mirna_id %in% AD$mirna_id]         <- "AD"
covs$RESIL[covs$mirna_id %in% normals$mirna_id]    <- "Control"
covs$RESIL[covs$mirna_id %in% presym$mirna_id]     <- "Presym"

covs$Diagnosis <- covs$RESIL




# Coerce variables of interest
# msex: 1 = male, 0 = female (ROSMAP coding). APOE4 carrier = any e4 allele.
covs <- covs %>%
  mutate(
    braaksc   = factor(braaksc, ordered = TRUE),
    ceradsc   = factor(ceradsc, ordered = TRUE),
    cogdx     = factor(cogdx,   ordered = TRUE),
    Sex       = factor(ifelse(msex == 1, "Male", "Female"), levels = c("Female", "Male")),
    APOE4     = factor(ifelse(apoe_genotype %in% c(24, 34, 44), "e4 carrier", "Non-carrier"),
                       levels = c("Non-carrier", "e4 carrier")),
    Diagnosis = factor(Diagnosis, levels = c("Control", "Presym", "Resilient", "AD", "OTHER"))
  )

table(covs$Diagnosis, useNA = "ifany")




#### 1. Cohort characteristics tables ####

summ_vars <- c("age_death", "Sex", "APOE4", "educ", "pmi", "RIN",
               "braaksc", "ceradsc", "cogdx")

summ_labels <- list(age_death = "Age at death (years)",
                    educ      = "Education (years)",
                    pmi       = "Post-mortem interval (h)",
                    RIN       = "RNA integrity number",
                    braaksc   = "Braak stage",
                    ceradsc   = "CERAD score",
                    cogdx     = "Clinical diagnosis (cogdx)")

make_tab <- function(df, by = "Diagnosis") {
  vars <- setdiff(summ_vars, by)
  tab <- df |>
    dplyr::select(all_of(c(by, vars))) |>
    droplevels() |>
    tbl_summary(
      by = all_of(by),
      type = list(c(age_death, educ, pmi, RIN) ~ "continuous2"),
      statistic = list(all_continuous2() ~ c("{mean} ± {sd}",
                                             "{median} [{p25}, {p75}]",
                                             "{min} – {max}"),
                       all_categorical() ~ "{n} ({p}%)"),
      digits = list(all_continuous2() ~ 1),
      label = summ_labels[names(summ_labels) %in% vars],
      missing_text = "Missing"
    ) |>
    add_overall()
  # add_p() needs the cardx package; group tests are also computed in base R below
  if (requireNamespace("cardx", quietly = TRUE)) {
    # pathology/cogdx define the groups, so only test demographics
    tab <- tab |>
      add_p(include = c(age_death, Sex, APOE4, educ, pmi, RIN),
            test = list(all_continuous() ~ "kruskal.test",
                        all_categorical() ~ "chisq.test"))
  }
  tab
}

# Demographic group differences (base R): across the four defined groups
# (Kruskal-Wallis / chi-square) and for the main contrast, AD vs Resilient
# (Wilcoxon / chi-square)
group_tests <- function(df, groups, label) {
  d <- droplevels(df[df$Diagnosis %in% groups, ])
  cont <- c("age_death", "educ", "pmi", "RIN")
  catg <- c("Sex", "APOE4")
  bind_rows(
    map_dfr(cont, ~ tibble(Variable = .x,
                           Test = if (length(groups) > 2) "Kruskal-Wallis" else "Wilcoxon",
                           P.Value = if (length(groups) > 2)
                             kruskal.test(d[[.x]] ~ d$Diagnosis)$p.value
                           else wilcox.test(d[[.x]] ~ d$Diagnosis)$p.value)),
    map_dfr(catg, ~ tibble(Variable = .x, Test = "Chi-square",
                           P.Value = suppressWarnings(chisq.test(table(d[[.x]], d$Diagnosis))$p.value)))
  ) %>%
    mutate(Groups = label, P.Value = signif(P.Value, 3), .before = 1)
}

save_tab <- function(tab, name) {
  write_csv(gtsummary::as_tibble(tab, col_labels = TRUE),
            file.path(meta.dir, paste0(name, ".csv")))
  if (requireNamespace("flextable", quietly = TRUE)) {
    flextable::save_as_docx(as_flex_table(tab),
                            path = file.path(meta.dir, paste0(name, ".docx")))
  }
  invisible(tab)
}

# a) All miRNA samples (n = 568), by group
tab_all <- make_tab(covs)
save_tab(tab_all, "Table_cohort_miRNA_all")
if (interactive()) tab_all

demo_tests <- bind_rows(
  group_tests(covs, c("Control", "Presym", "Resilient", "AD"), "Control/Presym/Resilient/AD"),
  group_tests(covs, c("Resilient", "AD"), "AD vs Resilient")
)
write_csv(demo_tests, file.path(meta.dir, "Demographics_group_tests.csv"))
print(as.data.frame(demo_tests))

# b) Samples with matched miRNA + mRNA
mRNA_EXPES <- readRDS("output/DEGenes_bulk/RESIDUAL_gene_expression.RDS")
covs_matched <- covs[covs$Sample_ID %in% colnames(mRNA_EXPES), ]
rm(mRNA_EXPES)

tab_matched <- make_tab(covs_matched)
save_tab(tab_matched, "Table_cohort_miRNA_mRNA_matched")
if (interactive()) tab_matched

# c) Proportional-odds subset (cogdx 1, 2, 4; see 03_01_miRNA_PathologyOddsModel.R)
covs_podds <- covs %>% filter(cogdx %in% c(1, 2, 4))
tab_podds  <- make_tab(covs_podds, by = "cogdx")
save_tab(tab_podds, "Table_cohort_propOdds_by_cogdx")

# d) Co-expression subset (matched miRNA/mRNA AND cogdx 1, 2, 4)
covs_coexp <- covs_matched %>% filter(cogdx %in% c(1, 2, 4))
tab_coexp  <- make_tab(covs_coexp, by = "cogdx")
save_tab(tab_coexp, "Table_cohort_coexpression_by_cogdx")




#### 2. Plain numbers for the Abstract (n, age, sex) ####

demo_line <- function(df, label) {
  tibble(
    Set          = label,
    N            = nrow(df),
    Age_mean     = mean(df$age_death, na.rm = TRUE),
    Age_sd       = sd(df$age_death, na.rm = TRUE),
    Age_median   = median(df$age_death, na.rm = TRUE),
    Age_IQR_low  = quantile(df$age_death, 0.25, na.rm = TRUE),
    Age_IQR_high = quantile(df$age_death, 0.75, na.rm = TRUE),
    Age_min      = min(df$age_death, na.rm = TRUE),
    Age_max      = max(df$age_death, na.rm = TRUE),
    Female_n     = sum(df$Sex == "Female", na.rm = TRUE),
    Female_pct   = 100 * mean(df$Sex == "Female", na.rm = TRUE),
    Male_n       = sum(df$Sex == "Male", na.rm = TRUE),
    Male_pct     = 100 * mean(df$Sex == "Male", na.rm = TRUE),
    APOE4_pct    = 100 * mean(df$APOE4 == "e4 carrier", na.rm = TRUE)
  )
}

demo_tab <- bind_rows(
  demo_line(covs, "All miRNA samples"),
  map_dfr(levels(droplevels(covs$Diagnosis)),
          ~ demo_line(covs[covs$Diagnosis == .x, ], paste0("miRNA: ", .x))),
  demo_line(covs_matched, "Matched miRNA + mRNA"),
  demo_line(covs_podds,   "Proportional-odds subset (cogdx 1/2/4)"),
  demo_line(covs_coexp,   "Co-expression subset (matched, cogdx 1/2/4)")
) %>%
  mutate(across(where(is.double), ~ round(.x, 1)))

write_csv(demo_tab, file.path(meta.dir, "Abstract_demographics.csv"))
print(as.data.frame(demo_tab))

# Sex breakdown per group (supports the sex-specific analysis, e.g. males: AD vs Resilient)
sex_by_group <- as.data.frame.matrix(table(covs$Diagnosis, covs$Sex))
sex_by_group
write.csv(sex_by_group, file.path(meta.dir, "Sex_by_group.csv"))




#### 3. Statistical threshold and effect estimates for the main findings ####

# Differential expression (limma, adjusted for sex, age, PMI, APOE, RIN; BH-FDR).
# Expression is log2, so logFC is log2 fold change. The 95% CI uses SE = logFC / t
# with a normal approximation (moderated t df is large, > 250).
summarize_de <- function(file, label) {
  path <- file.path(demir.dir, file)
  if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
  read_csv(path, show_col_types = FALSE, name_repair = "unique_quiet") %>%
    rename(miRNA = 1) %>%
    mutate(Comparison = label,
           SE         = logFC / t,
           CI_low     = logFC - qnorm(0.975) * SE,
           CI_high    = logFC + qnorm(0.975) * SE,
           FoldChange = 2^logFC)
}

de_files <- c(
  "AD vs Resilient"              = "deMirsADvsResilient.csv",
  "AD vs Control"                = "deMirsADvsControl.csv",
  "Resilient vs Control"         = "deMirsResilientvsNormal.csv",
  "AD vs Presym"                 = "deMirsADvsPresym.csv",
  "Resilient vs Presym"          = "deMirsResilientvsPresym.csv",
  "Presym vs Control"            = "deMirsPresymvsControl.csv",
  "AD vs Resilient (males)"      = "deMirs_sexSpec_AD_vs_Resilient_in_Males.csv",
  "AD vs Resilient (females)"    = "deMirs_sexSpec_AD_vs_Resilient_in_Females.csv",
  "AD vs Control (males)"        = "deMirs_sexSpec_AD_vs_Ctrl_in_Males.csv",
  "AD vs Control (females)"      = "deMirs_sexSpec_AD_vs_Ctrl_in_Females.csv"
)

de_all <- imap_dfr(de_files, ~ summarize_de(.x, .y))

# Number tested / significant per comparison, and the nominal p-value that
# corresponds to FDR < 0.05 (largest raw p among significant miRNAs)
de_overview <- de_all %>%
  group_by(Comparison) %>%
  summarise(
    n_tested       = n(),
    n_sig_FDR05    = sum(adj.P.Val < fdr.cut),
    n_up           = sum(adj.P.Val < fdr.cut & logFC > 0),
    n_down         = sum(adj.P.Val < fdr.cut & logFC < 0),
    nominal_p_at_FDR05 = if (any(adj.P.Val < fdr.cut)) max(P.Value[adj.P.Val < fdr.cut]) else NA_real_,
    min_FDR        = min(adj.P.Val),
    abs_log2FC_range_sig = if (any(adj.P.Val < fdr.cut))
      paste(signif(range(abs(logFC[adj.P.Val < fdr.cut])), 2), collapse = " – ") else NA_character_,
    .groups = "drop"
  ) %>%
  arrange(match(Comparison, names(de_files)))

write_csv(de_overview, file.path(meta.dir, "DE_overview_FDR05.csv"))
print(as.data.frame(de_overview))

de_sig <- de_all %>%
  filter(adj.P.Val < fdr.cut) %>%
  transmute(Comparison, miRNA,
            log2FC = round(logFC, 3),
            CI95   = sprintf("%.2f to %.2f", CI_low, CI_high),
            FoldChange = round(FoldChange, 2),
            P.Value = signif(P.Value, 3),
            FDR     = signif(adj.P.Val, 3))

write_csv(de_sig, file.path(meta.dir, "DE_significant_effect_sizes.csv"))
print(as.data.frame(de_sig))


# Proportional odds models (VGAM cumulative, parallel; 03_01_miRNA_PathologyOddsModel.R).
# VGAM models logit P(Y <= j), so OR for a MORE severe category per 1-unit (log2)
# increase in residual miRNA expression = exp(-Estimate). CERAD was reversed
# (4 - ceradsc) so higher = more plaque for all three outcomes.
summarize_podds <- function(file, label) {
  path <- file.path(demir.dir, file)
  if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
  read_csv(path, show_col_types = FALSE, name_repair = "unique_quiet") %>%
    transmute(Outcome = label,
              miRNA   = gene,
              Estimate,
              SE      = `Std. Error`,
              OR_more_severe = exp(-Estimate),
              OR_CI_low      = exp(-(Estimate + qnorm(0.975) * SE)),
              OR_CI_high     = exp(-(Estimate - qnorm(0.975) * SE)),
              P.Value = pval,
              FDR     = p.adj)
}

podds_files <- c(
  "Cognition (cogdx)"                     = "cognition_miRNA_prop_odds.csv",
  "Tangles (Braak)"                       = "tangle_miRNA_prop_odds.csv",
  "Plaques (CERAD)"                       = "plaque_miRNA_prop_odds.csv",
  "Cognition adjusted for plaques+tangles" = "adjusted_cognition_miRNA.csv"
)

podds_all <- imap_dfr(podds_files, ~ summarize_podds(.x, .y))

podds_overview <- podds_all %>%
  group_by(Outcome) %>%
  summarise(n_tested    = n(),
            n_sig_FDR05 = sum(FDR < fdr.cut),
            n_higher_risk = sum(FDR < fdr.cut & OR_more_severe > 1),
            n_lower_risk  = sum(FDR < fdr.cut & OR_more_severe < 1),
            .groups = "drop")
podds_overview <- podds_overview %>%
  mutate(n_subjects = nrow(covs_podds))
write_csv(podds_overview, file.path(meta.dir, "PropOdds_overview_FDR05.csv"))
print(as.data.frame(podds_overview))

podds_sig <- podds_all %>%
  filter(FDR < fdr.cut) %>%
  mutate(OR_CI95 = sprintf("%.2f (%.2f–%.2f)", OR_more_severe, OR_CI_low, OR_CI_high))
write_csv(podds_sig, file.path(meta.dir, "PropOdds_significant_effect_sizes.csv"))




#### 4. Key miRNAs named in the Abstract ####

key_mirs <- c("hsa-miR-132-3p", "hsa-miR-129-5p", "hsa-miR-362-3p", "hsa-miR-433-3p",
              "hsa-miR-335-5p", "hsa-miR-19b-3p", "hsa-miR-199a-5p", "hsa-miR-7-5p")

key_de <- de_all %>%
  filter(miRNA %in% key_mirs,
         Comparison %in% c("AD vs Resilient", "AD vs Control",
                           "AD vs Resilient (males)", "AD vs Resilient (females)")) %>%
  transmute(miRNA, Comparison,
            log2FC_CI = sprintf("%.2f (%.2f to %.2f)", logFC, CI_low, CI_high),
            FDR = signif(adj.P.Val, 2)) %>%
  arrange(miRNA, Comparison)

key_podds <- podds_all %>%
  filter(miRNA %in% key_mirs) %>%
  transmute(miRNA, Outcome,
            OR_CI = sprintf("%.2f (%.2f–%.2f)", OR_more_severe, OR_CI_low, OR_CI_high),
            FDR = signif(FDR, 2)) %>%
  arrange(miRNA, Outcome)

write_csv(key_de,    file.path(meta.dir, "KeyMiRNAs_DE.csv"))
write_csv(key_podds, file.path(meta.dir, "KeyMiRNAs_PropOdds.csv"))
print(as.data.frame(key_de))
print(as.data.frame(key_podds))




#### 5. Draft sentences for the Abstract ####

fmt_age <- function(df) {
  sprintf("%.1f ± %.1f years (range %.1f–%.1f)",
          mean(df$age_death, na.rm = TRUE), sd(df$age_death, na.rm = TRUE),
          min(df$age_death, na.rm = TRUE),  max(df$age_death, na.rm = TRUE))
}
fmt_sex <- function(df) {
  sprintf("%.1f%% female", 100 * mean(df$Sex == "Female", na.rm = TRUE))
}
grp <- function(g) covs[covs$Diagnosis == g, ]

ad_res <- de_overview %>% filter(Comparison == "AD vs Resilient")

abstract_txt <- c(
  sprintf("Total miRNA samples: n = %d; age at death %s; %s (%d females, %d males).",
          nrow(covs), fmt_age(covs), fmt_sex(covs),
          sum(covs$Sex == "Female"), sum(covs$Sex == "Male")),
  map_chr(c("Resilient", "AD", "Control", "Presym"), function(g) {
    d <- grp(g)
    sprintf("  %s: n = %d; age %s; %s.", g, nrow(d), fmt_age(d), fmt_sex(d))
  }),
  sprintf("Matched miRNA/mRNA: n = %d; age %s; %s.",
          nrow(covs_matched), fmt_age(covs_matched), fmt_sex(covs_matched)),
  sprintf("Significance threshold: Benjamini–Hochberg FDR < %.2f (AD vs Resilient: nominal p <= %.2g; %d of %d miRNAs significant, %d up / %d down in AD; |log2FC| %s).",
          fdr.cut, ad_res$nominal_p_at_FDR05, ad_res$n_sig_FDR05, ad_res$n_tested,
          ad_res$n_up, ad_res$n_down, ad_res$abs_log2FC_range_sig)
)

writeLines(abstract_txt, file.path(meta.dir, "Abstract_numbers.txt"))
cat(abstract_txt, sep = "\n")
