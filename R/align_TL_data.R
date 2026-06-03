#' Align target and source datasets
#'
#' Remove observations with missing outcomes, retain only common
#' features across all datasets, and align feature order.
#'
#' @param target A target data frame.
#' @param source A list of source data frames.
#' @param outcome Name of the outcome variable.
#'
#' @return A list containing:
#' \itemize{
#'   \item target: aligned target data frame.
#'   \item source: aligned source data frames.
#'   \item common_features: common retained features.
#'   \item removed_outcome_missing: removed observations due to missing outcome.
#' }
#'
#' @examples
#' \dontrun{
#' align <- align_TL_data(
#'   target = target_data,
#'   source = source_data,
#'   outcome = "Y"
#' )
#' }
#' @export
align_TL_data <- function(target, source, outcome) {

  if (!is.data.frame(target)) {
    stop("target must be a data.frame.")
  }

  if (!is.list(source) || !all(sapply(source, is.data.frame))) {
    stop("source must be a list of data.frames.")
  }

  all_data <- c(list(target = target), source)

  # Check outcome exists
  for (i in seq_along(all_data)) {
    if (!outcome %in% colnames(all_data[[i]])) {
      stop(paste0(
        "Outcome column '", outcome,
        "' is missing in dataset ", i, "."
      ))
    }
  }

  # Remove observations with missing outcomes
  removed_outcome_missing <- vector("list", length(all_data))

  for (i in seq_along(all_data)) {
    dat <- all_data[[i]]

    removed_outcome_missing[[i]] <- which(is.na(dat[[outcome]]))
    dat <- dat[!is.na(dat[[outcome]]), , drop = FALSE]

    all_data[[i]] <- dat
  }

  # Move outcome to first column
  all_data <- lapply(all_data, function(dat) {
    dat[, c(outcome, setdiff(colnames(dat), outcome)), drop = FALSE]
  })

  # Find common features excluding outcome
  feature_sets <- lapply(all_data, function(dat) {
    setdiff(colnames(dat), outcome)
  })

  common_features <- Reduce(intersect, feature_sets)

  # Use target feature order
  target_feature_order <- setdiff(
    colnames(all_data[[1]]),
    outcome
  )

  common_features <- target_feature_order[
    target_feature_order %in% common_features
  ]

  # Align all datasets
  all_data <- lapply(all_data, function(dat) {
    dat[, c(outcome, common_features), drop = FALSE]
  })

  target_aligned <- all_data[[1]]
  source_aligned <- all_data[-1]

  if (is.null(names(source)) || any(names(source) == "")) {
    names(source_aligned) <- paste0("source", seq_along(source_aligned))
  } else {
    names(source_aligned) <- names(source)
  }

  return(list(
    target = target_aligned,
    source = source_aligned,
    common_features = common_features,
    removed_outcome_missing = removed_outcome_missing
  ))
}
