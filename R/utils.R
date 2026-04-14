
.check_dataset <- function(x, name = deparse(substitute(x))) {
  if (is.null(x) || !is.list(x) || !all(c("x", "y") %in% names(x))) {
    stop(sprintf("`%s` must be a list with components `x` and `y`.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.check_source_list <- function(x, name = deparse(substitute(x))) {
  if (is.null(x) || !is.list(x) || length(x) == 0) {
    stop(sprintf("`%s` must be a non-empty list of datasets.", name), call. = FALSE)
  }
  for (i in seq_along(x)) .check_dataset(x[[i]], sprintf("%s[[%d]]", name, i))
  invisible(TRUE)
}

.check_source_id <- function(source_id) {
  if (!source_id %in% c("auto", "all")) {
    stop('`source_id` must be either "auto" or "all".')
  }
  invisible(TRUE)
}

.check_beta_hat <- function(beta.hat) {
  if (is.null(beta.hat)) {
    stop("`beta.hat` must be supplied.", call. = FALSE)
  }
  invisible(TRUE)
}
