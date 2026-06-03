#' Process microbiome transfer learning datasets
#'
#' Process aligned microbiome datasets by removing highly-missing samples
#' and features, converting microbiome variables to log-relative abundance,
#' and imputing non-microbiome covariates.
#'
#' @param target Aligned target data frame.
#' @param source List of aligned source data frames.
#' @param outcome Outcome variable name.
#' @param covariate Vector of non-microbiome covariate names.
#' @param microbiome.scale One of \code{"count"}, \code{"relative_abundance"},
#'   or \code{"log_relative_abundance"}.
#' @param microbe_missing_rate Missingness threshold for microbiomes.
#' @param covariate_missing_rate Missingness threshold for covariates.
#' @param sample_missing_rate Missingness threshold for observations.
#' @param pseudo.count Small positive value replacing zeros. If \code{NULL},
#'   half of the minimum positive relative abundance is used, bounded below by
#'   \code{1e-10}.
#' @param use.missForest Logical indicating whether missForest is used.
#'
#' @return A list containing processed target and source datasets.
#'
#' @examples
#' \dontrun{
#' processed <- process_microbiome_TL_data(
#'   target = target_aligned,
#'   source = source_aligned,
#'   outcome = "Y",
#'   covariate = c("cov1", "cov2")
#' )
#' }
#'
#' @export
process_microbiome_TL_data <- function(
    target,
    source,
    outcome,
    covariate,
    microbiome.scale = c("count", "relative_abundance", "log_relative_abundance"),
    microbe_missing_rate = 0.2,
    covariate_missing_rate = 0.5,
    sample_missing_rate = 0.5,
    pseudo.count = NULL,
    use.missForest = TRUE
) {

  microbiome.scale <- match.arg(microbiome.scale)

  if (use.missForest && !requireNamespace("missForest", quietly = TRUE)) {
    stop("Please install missForest first: install.packages('missForest')")
  }

  if (!is.data.frame(target)) {
    stop("target must be a data.frame.")
  }

  if (!is.list(source) || !all(sapply(source, is.data.frame))) {
    stop("source must be a list of data.frames.")
  }

  all_data <- c(list(target = target), source)

  for (k in seq_along(all_data)) {
    if (!outcome %in% colnames(all_data[[k]])) {
      stop(paste0("Outcome column is missing in dataset ", k, "."))
    }
  }

  processed_each <- vector("list", length(all_data))
  dropped_samples <- vector("list", length(all_data))
  dropped_covariates <- vector("list", length(all_data))
  dropped_microbes <- vector("list", length(all_data))

  for (k in seq_along(all_data)) {

    dat <- all_data[[k]]

    cov_k <- intersect(covariate, colnames(dat))
    microbe_k <- setdiff(colnames(dat), c(outcome, covariate))

    dat <- dat[, c(outcome, cov_k, microbe_k), drop = FALSE]

    predictor_names <- c(cov_k, microbe_k)

    if (length(predictor_names) == 0) {
      stop(paste0("No predictors are available in dataset ", k, "."))
    }

    miss_mat <- is.na(dat[, predictor_names, drop = FALSE])
    sample_miss <- rowMeans(miss_mat)
    keep_sample <- sample_miss <= sample_missing_rate

    dropped_samples[[k]] <- which(!keep_sample)
    dat <- dat[keep_sample, , drop = FALSE]

    if (nrow(dat) == 0) {
      stop(paste0("Dataset ", k, " has 0 rows after sample filtering."))
    }

    if (length(cov_k) > 0) {
      cov_miss <- colMeans(is.na(dat[, cov_k, drop = FALSE]))
      keep_cov <- cov_miss <= covariate_missing_rate
      dropped_covariates[[k]] <- names(cov_miss)[!keep_cov]
      cov_k <- names(cov_miss)[keep_cov]
    } else {
      dropped_covariates[[k]] <- character(0)
    }

    if (length(microbe_k) > 0) {
      microbe_miss <- colMeans(is.na(dat[, microbe_k, drop = FALSE]))
      keep_microbe <- microbe_miss <= microbe_missing_rate
      dropped_microbes[[k]] <- names(microbe_miss)[!keep_microbe]
      microbe_k <- names(microbe_miss)[keep_microbe]
    } else {
      dropped_microbes[[k]] <- character(0)
    }

    processed_each[[k]] <- dat[, c(outcome, cov_k, microbe_k), drop = FALSE]
  }

  cov_sets <- lapply(processed_each, function(dat) {
    intersect(covariate, colnames(dat))
  })

  microbe_sets <- lapply(processed_each, function(dat) {
    setdiff(colnames(dat), c(outcome, covariate))
  })

  common_cov <- Reduce(intersect, cov_sets)
  common_microbe <- Reduce(intersect, microbe_sets)

  common_cov <- covariate[covariate %in% common_cov]

  target_microbe_order <- setdiff(
    colnames(processed_each[[1]]),
    c(outcome, covariate)
  )

  common_microbe <- target_microbe_order[
    target_microbe_order %in% common_microbe
  ]

  if (length(common_microbe) == 0) {
    stop("No common microbiome features remain after filtering.")
  }

  processed_each <- lapply(processed_each, function(dat) {
    dat[, c(outcome, common_cov, common_microbe), drop = FALSE]
  })

  dataset_id <- rep(seq_along(processed_each), sapply(processed_each, nrow))
  combined <- do.call(rbind, processed_each)

  Y <- combined[[outcome]]
  cov_dat <- combined[, common_cov, drop = FALSE]
  microb_dat <- as.matrix(combined[, common_microbe, drop = FALSE])
  storage.mode(microb_dat) <- "numeric"

  if (microbiome.scale == "count") {

    microb_dat[is.na(microb_dat)] <- 0
    microb_dat[microb_dat < 0] <- 0

    row_sum <- rowSums(microb_dat)
    keep_nonzero <- row_sum > 0

    if (!any(keep_nonzero)) {
      stop("All samples have zero microbiome counts after preprocessing.")
    }

    Y <- Y[keep_nonzero]
    cov_dat <- cov_dat[keep_nonzero, , drop = FALSE]
    microb_dat <- microb_dat[keep_nonzero, , drop = FALSE]
    dataset_id <- dataset_id[keep_nonzero]

    RA <- microb_dat / rowSums(microb_dat)
  }

  if (microbiome.scale == "relative_abundance") {

    microb_dat[is.na(microb_dat)] <- 0
    microb_dat[microb_dat < 0] <- 0

    row_sum <- rowSums(microb_dat)
    keep_nonzero <- row_sum > 0

    if (!any(keep_nonzero)) {
      stop("All samples have zero microbiome relative abundance after preprocessing.")
    }

    Y <- Y[keep_nonzero]
    cov_dat <- cov_dat[keep_nonzero, , drop = FALSE]
    microb_dat <- microb_dat[keep_nonzero, , drop = FALSE]
    dataset_id <- dataset_id[keep_nonzero]

    RA <- microb_dat / rowSums(microb_dat)
  }

  if (microbiome.scale == "log_relative_abundance") {

    RA_raw <- exp(microb_dat)
    RA_raw[is.na(RA_raw)] <- 0
    RA_raw[RA_raw < 0] <- 0

    row_sum <- rowSums(RA_raw)
    keep_nonzero <- row_sum > 0

    if (!any(keep_nonzero)) {
      stop("All samples have zero microbiome relative abundance after converting from log scale.")
    }

    Y <- Y[keep_nonzero]
    cov_dat <- cov_dat[keep_nonzero, , drop = FALSE]
    RA_raw <- RA_raw[keep_nonzero, , drop = FALSE]
    dataset_id <- dataset_id[keep_nonzero]

    RA <- RA_raw / rowSums(RA_raw)
  }

  if (is.null(pseudo.count)) {
    positive_RA <- RA[RA > 0]

    if (length(positive_RA) == 0) {
      stop("No positive relative abundance values are available for pseudo-count selection.")
    }

    pseudo.count <- max(
      min(positive_RA, na.rm = TRUE) / 2,
      1e-10
    )
  }

  RA[RA == 0] <- pseudo.count
  RA <- RA / rowSums(RA)

  logRA_df <- as.data.frame(log(RA), check.names = FALSE)
  colnames(logRA_df) <- common_microbe

  cov_imp <- cov_dat

  if (use.missForest && ncol(cov_imp) > 0 && any(is.na(cov_imp))) {

    X_for_impute <- data.frame(
      cov_imp,
      logRA_df,
      check.names = FALSE
    )

    for (v in common_cov) {
      obs_vals <- stats::na.omit(unique(X_for_impute[[v]]))

      if (length(obs_vals) <= 2 && all(obs_vals %in% c(0, 1))) {
        X_for_impute[[v]] <- factor(X_for_impute[[v]], levels = c(0, 1))
      }
    }

    X_imp <- missForest::missForest(X_for_impute)$ximp
    cov_imp <- X_imp[, common_cov, drop = FALSE]
  }

  combined_processed <- data.frame(
    cov_imp,
    logRA_df,
    check.names = FALSE
  )

  tmp_df <- data.frame(
    dataset_id = dataset_id,
    outcome_tmp = Y,
    combined_processed,
    check.names = FALSE
  )

  colnames(tmp_df)[2] <- outcome

  processed_split <- split(tmp_df, tmp_df$dataset_id)

  target_df <- processed_split[[as.character(1)]]

  target_data <- list(
    y = target_df[[outcome]],
    x = target_df[, setdiff(colnames(target_df), c("dataset_id", outcome)), drop = FALSE]
  )

  source_data <- vector("list", length(source))

  for (k in seq_along(source)) {

    source_df <- processed_split[[as.character(k + 1)]]

    if (is.null(source_df)) {
      stop(paste0("Source dataset ", k, " has no samples after preprocessing."))
    }

    source_data[[k]] <- list(
      y = source_df[[outcome]],
      x = source_df[, setdiff(colnames(source_df), c("dataset_id", outcome)), drop = FALSE]
    )
  }

  if (is.null(names(source)) || any(names(source) == "")) {
    names(source_data) <- paste0("source", seq_along(source_data))
  } else {
    names(source_data) <- names(source)
  }

  return(list(
    target_data = target_data,
    source_data = source_data,
    outcome = outcome,
    covariates = common_cov,
    microbiome_features = common_microbe,
    dropped_samples = dropped_samples,
    dropped_covariates = dropped_covariates,
    dropped_microbes = dropped_microbes,
    pseudo.count = pseudo.count
  ))
}
