#!/usr/bin/env Rscript

# Cox proportional-hazards analyses for organ risk indices.
#
# Default inputs in the working directory:
#   organ_risk_index.csv and sample_metadata.csv
# Optional usage:
#   Rscript cox-analysis.R [ORGAN_INDEX.csv] [SAMPLE_METADATA.csv] [OUTPUT_DIR]

suppressPackageStartupMessages({
  library(survival)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3) {
  stop("Usage: Rscript cox-analysis.R [ORGAN_INDEX.csv] [SAMPLE_METADATA.csv] [OUTPUT_DIR]")
}

organ_file <- if (length(args) >= 1) args[[1]] else "organ_risk_index.csv"
metadata_file <- if (length(args) >= 2) args[[2]] else "sample_metadata.csv"
output_dir <- if (length(args) == 3) args[[3]] else "cox-results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

organ_raw <- read.csv(organ_file, row.names = 1, check.names = FALSE)
info_raw <- read.csv(metadata_file, row.names = 1, check.names = FALSE)

required <- c("Group", "Year")
missing_required <- setdiff(required, colnames(info_raw))
if (length(missing_required) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_required, collapse = ", "))
}

common_ids <- intersect(rownames(org_raw), rownames(info_raw))
if (length(common_ids) < 3) {
  stop("At least three shared sample IDs are required.")
}
org_raw <- org_raw[common_ids, , drop = FALSE]
info <- info_raw[common_ids, , drop = FALSE]
org_raw[] <- lapply(org_raw, function(x) as.numeric(as.character(x)))
org_raw <- org_raw[, vapply(org_raw, function(x) any(is.finite(x)), logical(1)), drop = FALSE]
if (ncol(org_raw) == 0) stop("No numeric organ-risk columns were found.")
colnames(org_raw) <- make.names(colnames(org_raw), unique = TRUE)

info$Group <- factor(info$Group, levels = c("NC", "AD"))
if (anyNA(info$Group)) stop("Group must contain only 'NC' and 'AD'.")
info$E4_status <- if ("E4_status" %in% names(info)) {
  factor(info$E4_status, levels = c("E4_nocarrier", "E4_carrier"))
} else NULL
if ("Sex" %in% names(info)) info$Sex <- factor(info$Sex)
if ("HQ" %in% names(info)) info$HQ <- factor(info$HQ)
if ("Range" %in% names(info)) info$Range <- factor(info$Range)

info$time <- pmax(as.numeric(info$Year), 1e-3)
info$event <- as.integer(info$Group == "AD")
org <- as.data.frame(scale(org_raw))

safe_cox <- function(formula, data) {
  tryCatch(
    suppressWarnings(coxph(formula, data = data, ties = "efron", x = TRUE)),
    error = function(e) NULL
  )
}

empty_main_row <- function(organ, n, events) {
  data.frame(organ = organ, beta = NA_real_, se = NA_real_, z = NA_real_,
             p = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
             n = n, nevent = events)
}

main_covars <- intersect(c("E4_status", "Age", "Sex", "BMI", "TDI", "HQ"), names(info))
run_main_per_organ <- function() {
  rows <- lapply(names(org), function(organ_name) {
    dat <- cbind(info, organ_value = org[[organ_name]])
    rhs <- c("organ_value", main_covars)
    dat <- dat[complete.cases(dat[, c("time", "event", rhs), drop = FALSE]), , drop = FALSE]
    fit <- safe_cox(as.formula(paste("Surv(time, event) ~", paste(rhs, collapse = " + "))), dat)
    if (is.null(fit) || !"organ_value" %in% names(coef(fit))) {
      return(empty_main_row(organ_name, nrow(dat), sum(dat$event)))
    }
    sm <- summary(fit)$coef["organ_value", , drop = FALSE]
    ci <- tryCatch(confint(fit)["organ_value", ], error = function(e) c(NA_real_, NA_real_))
    data.frame(organ = organ_name, beta = sm[1, "coef"], se = sm[1, "se(coef)"],
               z = sm[1, "z"], p = sm[1, "Pr(>|z|)"], HR = exp(sm[1, "coef"]),
               CI_low = exp(ci[1]), CI_high = exp(ci[2]), n = nrow(dat),
               nevent = sum(dat$event))
  })
  result <- do.call(rbind, rows)
  result$p_adj <- p.adjust(result$p, method = "BH")
  result[order(result$p_adj, result$p, na.last = TRUE), , drop = FALSE]
}

write.csv(run_main_per_organ(), file.path(output_dir, "cox_main_per_organ.csv"), row.names = FALSE)

if ("E4_status" %in% names(info) && nlevels(droplevels(info$E4_status)) > 1) {
  run_e4_interaction <- function(organ_name) {
    dat <- cbind(info, organ_value = org[[organ_name]])
    covars <- setdiff(main_covars, "E4_status")
    rhs <- c("organ_value * E4_status", covars)
    vars <- c("time", "event", "organ_value", "E4_status", covars)
    dat <- dat[complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE]
    fit <- safe_cox(as.formula(paste("Surv(time, event) ~", paste(rhs, collapse = " + "))), dat)
    coefficients <- if (!is.null(fit)) coef(fit) else NULL
    vc <- if (!is.null(fit)) tryCatch(vcov(fit), error = function(e) NULL) else NULL
    interaction_term <- grep("^organ_value:E4_status", names(coefficients), value = TRUE)[1]
    if (is.null(coefficients) || is.na(interaction_term) || is.null(vc) ||
        !"organ_value" %in% names(coefficients)) {
      return(data.frame(organ = organ_name, HR_noncarrier = NA, HR_carrier = NA,
                        HR_ratio_carrier_vs_noncarrier = NA, p_ratio = NA))
    }
    beta0 <- coefficients[["organ_value"]]
    se0 <- sqrt(vc["organ_value", "organ_value"])
    beta_int <- coefficients[[interaction_term]]
    se_int <- sqrt(vc[interaction_term, interaction_term])
    beta_carrier <- beta0 + beta_int
    se_carrier <- sqrt(vc["organ_value", "organ_value"] + vc[interaction_term, interaction_term] +
                       2 * vc["organ_value", interaction_term])
    data.frame(
      organ = organ_name,
      HR_noncarrier = exp(beta0),
      CI_low_noncarrier = exp(beta0 + c(-1, 1)[1] * 1.96 * se0),
      CI_high_noncarrier = exp(beta0 + c(-1, 1)[2] * 1.96 * se0),
      HR_carrier = exp(beta_carrier),
      CI_low_carrier = exp(beta_carrier - 1.96 * se_carrier),
      CI_high_carrier = exp(beta_carrier + 1.96 * se_carrier),
      HR_ratio_carrier_vs_noncarrier = exp(beta_int),
      p_ratio = 2 * pnorm(-abs(beta_int / se_int))
    )
  }
  e4_result <- do.call(rbind, lapply(names(org), run_e4_interaction))
  e4_result$p_adj_ratio <- p.adjust(e4_result$p_ratio, method = "BH")
  write.csv(e4_result, file.path(output_dir, "cox_e4_interaction.csv"), row.names = FALSE)
}

if ("Range" %in% names(info) && nlevels(droplevels(info$Range)) > 1) {
  stage_rows <- lapply(names(org), function(organ_name) {
    dat <- cbind(info, organ_value = org[[organ_name]])
    rhs <- c("organ_value * Range", main_covars)
    vars <- unique(c("time", "event", "organ_value", "Range", main_covars))
    dat <- dat[complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE]
    fit <- safe_cox(as.formula(paste("Surv(time, event) ~", paste(rhs, collapse = " + "))), dat)
    if (is.null(fit)) return(NULL)
    b <- coef(fit)
    v <- tryCatch(vcov(fit), error = function(e) NULL)
    interaction_terms <- grep("^organ_value:Range", names(b), value = TRUE)
    if (is.null(v) || length(interaction_terms) == 0) return(NULL)
    data.frame(organ = organ_name, term = interaction_terms,
               beta = unname(b[interaction_terms]),
               se = sqrt(diag(v)[interaction_terms]))
  })
  stage_result <- do.call(rbind, stage_rows)
  if (!is.null(stage_result) && nrow(stage_result) > 0) {
    stage_result$z <- stage_result$beta / stage_result$se
    stage_result$p <- 2 * pnorm(-abs(stage_result$z))
    stage_result$p_adj <- ave(stage_result$p, sub("^.*:Range", "", stage_result$term),
                              FUN = function(x) p.adjust(x, method = "BH"))
    write.csv(stage_result, file.path(output_dir, "cox_range_interactions.csv"), row.names = FALSE)
  }
}

message("Cox analyses completed. Results written to: ", normalizePath(output_dir, mustWork = FALSE))
