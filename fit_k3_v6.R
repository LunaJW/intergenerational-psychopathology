rm(list = ls())
set.seed(123)

# ---- libraries ----
library(slca)
library(tidyverse)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)

# ---- paths ----
in_file <- "C:/Users/15280/OneDrive/文档/2025 Summer Project/data/mh_p_cbcl_tscores.xlsx"
out_dir <- "C:/Users/15280/OneDrive/文档/2025 Summer Project/outputs"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_fit <- file.path(out_dir, "fit_k3.rds")
out_data <- file.path(out_dir, "cbcl_wide_k3.rds")
out_info <- file.path(out_dir, "sessionInfo_fit_k3.txt")

# ---- data import ----
cbcl_long_raw <- read_xlsx(in_file)

# Raw sample size before any filtering
n_raw_subjects <- cbcl_long_raw %>%
  dplyr::distinct(participant_id) %>%
  nrow()

# ---- data preparation ----
cbcl_long <- cbcl_long_raw %>%
  dplyr::rename(
    subjectkey = participant_id,
    eventname  = session_id,
    cbcl_scr_syn_anxdep_t = mh_p_cbcl__synd__anxdep_tscore,
    cbcl_scr_syn_withdep_t = mh_p_cbcl__synd__wthdep_tscore,
    cbcl_scr_syn_somatic_t = mh_p_cbcl__synd__som_tscore,
    cbcl_scr_syn_social_t = mh_p_cbcl__synd__soc_tscore,
    cbcl_scr_syn_thought_t = mh_p_cbcl__synd__tho_tscore,
    cbcl_scr_syn_attention_t = mh_p_cbcl__synd__attn_tscore,
    cbcl_scr_syn_rulebreak_t = mh_p_cbcl__synd__rule_tscore,
    cbcl_scr_syn_aggressive_t = mh_p_cbcl__synd__aggr_tscore
  ) %>%
  dplyr::mutate(
    time = dplyr::case_when(
      eventname == "ses-00A" ~ "baseline",
      eventname == "ses-01A" ~ "1y",
      eventname == "ses-02A" ~ "2y",
      eventname == "ses-03A" ~ "3y",
      eventname == "ses-04A" ~ "4y",
      eventname == "ses-05A" ~ "5y",
      eventname == "ses-06A" ~ "6y",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(time)) %>%
  # Exclude 6-year follow-up to match the manuscript/Rmd analysis
  dplyr::filter(time != "6y") %>%
  dplyr::mutate(
    time = factor(
      time,
      levels = c("baseline", "1y", "2y", "3y", "4y", "5y")
    )
  )

n_after_timefilter_subjects <- cbcl_long %>%
  dplyr::distinct(subjectkey) %>%
  nrow()

cat("Raw unique participants:", n_raw_subjects, "\n")
cat("Participants after restricting to baseline-5y:", n_after_timefilter_subjects, "\n")
cat("Time suffixes included:", paste(levels(cbcl_long$time), collapse = ", "), "\n")
cat("Total observations after time filtering:", nrow(cbcl_long), "\n")

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

# ---- dichotomization ----
cbcl_long_bin <- cbcl_long %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(symptom_vars_t),
      ~ dplyr::if_else(.x >= 65, 1L, 0L, missing = NA_integer_)
    )
  )

# ---- wide format ----
cbcl_wide_bin <- cbcl_long_bin %>%
  dplyr::select(subjectkey, time, dplyr::all_of(symptom_vars_t)) %>%
  tidyr::pivot_wider(
    names_from  = time,
    values_from = dplyr::all_of(symptom_vars_t),
    names_glue  = "{.value}_{time}"
  )

n_wide_subjects <- nrow(cbcl_wide_bin)

# ---- filtering analytic sample ----
cbcl_wide_k3 <- cbcl_wide_bin %>%
  dplyr::mutate(
    non_missing_count = rowSums(
      !is.na(dplyr::across(dplyr::all_of(names(cbcl_wide_bin)[-1])))
    )
  ) %>%
  dplyr::filter(non_missing_count >= 3L * length(symptom_vars_t)) %>%
  dplyr::select(-non_missing_count) %>%
  dplyr::mutate(rowid = dplyr::row_number())

n_final_subjects <- nrow(cbcl_wide_k3)

cat("Wide-format participants before filtering:", n_wide_subjects, "\n")
cat("Final analytic sample:", n_final_subjects, "\n")

# ---- reorder wide columns by fixed time order ----
time_levels <- levels(cbcl_long$time)

ordered_symptom_cols <- unlist(
  lapply(time_levels, function(t) {
    paste0(symptom_vars_t, "_", t)
  }),
  use.names = FALSE
)

missing_ordered_cols <- setdiff(ordered_symptom_cols, names(cbcl_wide_k3))

if (length(missing_ordered_cols) > 0) {
  stop(
    "Missing expected ordered columns: ",
    paste(missing_ordered_cols, collapse = ", ")
  )
}

cbcl_wide_k3 <- cbcl_wide_k3 %>%
  dplyr::select(subjectkey, rowid, dplyr::all_of(ordered_symptom_cols))

# ---- save analysis dataset ----
saveRDS(cbcl_wide_k3, out_data)

# ---- sample flow table ----
sample_flow_table <- tibble::tibble(
  Step = c(
    "Raw unique participants in original file",
    "Participants after restricting to baseline-5y waves",
    "Wide-format participants before non-missing data filtering",
    "Final analytic sample with sufficient non-missing symptom data"
  ),
  N = c(
    n_raw_subjects,
    n_after_timefilter_subjects,
    n_wide_subjects,
    n_final_subjects
  )
)

print(sample_flow_table)

# ---- build formulas ----
time_suffixes <- levels(cbcl_long$time)

expected_cols <- unlist(
  lapply(time_suffixes, function(t) {
    paste0(symptom_vars_t, "_", t)
  }),
  use.names = FALSE
)

missing_cols <- setdiff(expected_cols, names(cbcl_wide_k3))

if (length(missing_cols) > 0) {
  stop(
    "Missing expected wide columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

L_names <- paste0("L", seq_along(time_suffixes))
K <- 3

cat("Time suffixes fixed order:", paste(time_suffixes, collapse = ", "), "\n")
cat("Latent variables:", paste(L_names, collapse = ", "), "\n")
cat("Number of classes:", K, "\n")

make_lv_formula <- function(L, k, time_suffix, symptoms) {
  rhs <- paste0(symptoms, "_", time_suffix, collapse = " + ")
  as.formula(paste0(L, "[", k, "] ~ ", rhs))
}

meas_forms_k3 <- Map(
  function(L, t) {
    make_lv_formula(L, K, t, symptom_vars_t)
  },
  L_names,
  time_suffixes
)

trans_forms <- lapply(
  seq_len(length(L_names) - 1),
  function(i) as.formula(paste0(L_names[i], " ~ ", L_names[i + 1]))
)

formulas_k3 <- c(meas_forms_k3, trans_forms)

cat("\nMeasurement formulas:\n")
print(meas_forms_k3)

cat("\nTransition formulas:\n")
print(trans_forms)

# ---- build model ----
model_k3 <- slca(
  formula = formulas_k3,
  data    = cbcl_wide_k3,
  control = slcaControl(
    structure = paste(L_names[-length(L_names)], "->", L_names[-1])
  )
)

message("Model built")
t1 <- Sys.time()

# ---- fit model ----
fit_k3 <- estimate(
  model_k3,
  data = cbcl_wide_k3,
  method = "em"
)

t2 <- Sys.time()
message("Estimate finished")
cat("Elapsed time:", as.character(t2 - t1), "\n")

# ---- save final objects ----
saveRDS(fit_k3, out_fit)
saveRDS(cbcl_wide_k3, out_data)

writeLines(
  capture.output(sessionInfo()),
  out_info
)