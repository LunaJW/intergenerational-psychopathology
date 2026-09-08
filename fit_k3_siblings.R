rm(list = ls())
set.seed(123)

# ---- libraries ----

library(slca)
library(tidyverse)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# ---- paths ----

in_data <- "C:/Users/15280/OneDrive/文档/2025 Summer Project/outputs/cbcl_wide_k3.rds"
dat_path <- "C:/Users/15280/OneDrive/文档/2025 Summer Project/latent transition analysis/ABCD_tables"
tsv_path <- file.path(dat_path, "ab_g_stc.tsv")

out_dir <- "C:/Users/15280/OneDrive/文档/2025 Summer Project/outputs"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_fit <- file.path(out_dir, "fit_k3_sibling.rds")
out_data <- file.path(out_dir, "cbcl_wide_k3_sibling.rds")
out_info <- file.path(out_dir, "sessionInfo_fit_k3_sibling.txt")

# ---- load original Aim 1 analytic dataset ----

cbcl_wide_k3 <- readRDS(in_data)

cat("Original Aim 1 N:", nrow(cbcl_wide_k3), "\n")

# ---- family information ----

stc_dat <- readr::read_tsv(
  tsv_path,
  show_col_types = FALSE
)

sibling_df <- stc_dat %>%
  dplyr::select(
    participant_id,
    ab_g_stc__design_id__fam__gen,
    ab_g_stc__design_famrel
  ) %>%
  dplyr::rename(subjectkey = participant_id) %>%
  dplyr::mutate(
    ab_g_stc__design_id__fam__gen =
      dplyr::na_if(ab_g_stc__design_id__fam__gen, "n/a")
  ) %>%
  dplyr::filter(!is.na(subjectkey)) %>%
  dplyr::distinct()

# ---- create one-participant-per-family sample ----

cbcl_wide_k3_sibling <- cbcl_wide_k3 %>%
  dplyr::left_join(
    sibling_df,
    by = "subjectkey"
  ) %>%
  dplyr::mutate(
    family_id_analysis = dplyr::if_else(
      is.na(ab_g_stc__design_id__fam__gen),
      subjectkey,
      ab_g_stc__design_id__fam__gen
    )
  )

set.seed(2025)

cbcl_wide_k3_sibling <- cbcl_wide_k3_sibling %>%
  dplyr::group_by(family_id_analysis) %>%
  dplyr::slice_sample(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    -ab_g_stc__design_id__fam__gen,
    -ab_g_stc__design_famrel,
    -family_id_analysis
  )

cat(
  "One-participant-per-family Aim 1 N:",
  nrow(cbcl_wide_k3_sibling),
  "\n"
)

cat(
  "Participants excluded:",
  nrow(cbcl_wide_k3) - nrow(cbcl_wide_k3_sibling),
  "\n"
)

# ---- symptom variables ----

symptom_vars_t <- c(
  "cbcl_scr_syn_anxdep_t",
  "cbcl_scr_syn_withdep_t",
  "cbcl_scr_syn_somatic_t",
  "cbcl_scr_syn_social_t",
  "cbcl_scr_syn_thought_t",
  "cbcl_scr_syn_attention_t",
  "cbcl_scr_syn_rulebreak_t",
  "cbcl_scr_syn_aggressive_t"
)

# ---- fixed time order ----

time_suffixes <- c(
  "baseline",
  "1y",
  "2y",
  "3y",
  "4y",
  "5y"
)

# ---- check expected columns ----

expected_cols <- unlist(
  lapply(time_suffixes, function(t) {
    paste0(symptom_vars_t, "_", t)
  }),
  use.names = FALSE
)

missing_cols <- setdiff(
  expected_cols,
  names(cbcl_wide_k3_sibling)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing expected wide columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ---- build formulas ----

L_names <- paste0(
  "L",
  seq_along(time_suffixes)
)

K <- 3

cat(
  "Time suffixes fixed order:",
  paste(time_suffixes, collapse = ", "),
  "\n"
)

cat(
  "Latent variables:",
  paste(L_names, collapse = ", "),
  "\n"
)

cat(
  "Number of classes:",
  K,
  "\n"
)

make_lv_formula <- function(L, k, time_suffix, symptoms) {
  rhs <- paste0(
    symptoms,
    "_",
    time_suffix,
    collapse = " + "
  )
  
  as.formula(
    paste0(
      L,
      "[",
      k,
      "] ~ ",
      rhs
    )
  )
}

meas_forms_k3 <- Map(
  function(L, t) {
    make_lv_formula(
      L,
      K,
      t,
      symptom_vars_t
    )
  },
  L_names,
  time_suffixes
)

trans_forms <- lapply(
  seq_len(length(L_names) - 1),
  function(i) {
    as.formula(
      paste0(
        L_names[i],
        " ~ ",
        L_names[i + 1]
      )
    )
  }
)

formulas_k3 <- c(
  meas_forms_k3,
  trans_forms
)

cat("\nMeasurement formulas:\n")
print(meas_forms_k3)

cat("\nTransition formulas:\n")
print(trans_forms)

# ---- build model ----

model_k3_sibling <- slca(
  formula = formulas_k3,
  data = cbcl_wide_k3_sibling,
  control = slcaControl(
    structure = paste(
      L_names[-length(L_names)],
      "->",
      L_names[-1]
    )
  )
)

message("Sibling sensitivity model built")

t1 <- Sys.time()

# ---- fit model ----

fit_k3_sibling <- estimate(
  model_k3_sibling,
  data = cbcl_wide_k3_sibling,
  method = "em"
)

t2 <- Sys.time()

message("Sibling sensitivity estimate finished")

cat(
  "Elapsed time:",
  as.character(t2 - t1),
  "\n"
)

# ---- save final objects ----

saveRDS(
  fit_k3_sibling,
  out_fit
)

saveRDS(
  cbcl_wide_k3_sibling,
  out_data
)

writeLines(
  capture.output(sessionInfo()),
  out_info
)

cat("\nSaved files:\n")
cat(out_fit, "\n")
cat(out_data, "\n")
cat(out_info, "\n")