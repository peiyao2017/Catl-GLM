utils::globalVariables(c(
  "beta_total_binomial",
  "beta_total_linear"
))
utils::globalVariables(c(
  "G_t_cpp",
  "G_t_debias_cpp",
  "G_t_ls_constrained_cpp",
  "G_t_ls_debias_constrained_cpp",
  "g_cpp",
  "g_ls_cpp",
  "grad_g_cpp",
  "grad_g_ls_cpp",
  "projker_ct_cpp",
  "proj_ker_ct_cpp",
  "prox_l1_shift_zerosum_rank1_cpp",
  "prox_l1_zerosum_rank1_cpp",
  "soft_thresh_cpp",
  "soft_thresh_debias_cpp"
))

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
