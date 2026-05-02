## ================================================================
##  Logistic transfer learning code (cleaned, NA-safe)
## ================================================================

info_logistic_n <- function(y, X, beta) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  beta <- as.numeric(beta)

  n <- length(y)

  ## beta always contains an intercept slot in position 1
  if (ncol(X) == length(beta) - 1) {
    X <- cbind(1, X)
  }

  if (nrow(X) != n) {
    stop("X and y dimension mismatch")
  }
  if (ncol(X) != length(beta)) {
    stop(sprintf(
      "beta length mismatch: ncol(X) = %d, length(beta) = %d",
      ncol(X), length(beta)
    ))
  }

  eta <- as.vector(X %*% beta)

  ## numerically stable logistic probability
  p <- stats::plogis(eta)
  w <- p * (1 - p)
  w[!is.finite(w)] <- 0

  WX <- X * w
  I <- t(X) %*% WX / n

  I
}

psi_logit <- function(u) {
  out <- numeric(length(u))
  large <- u > 0
  out[large] <- u[large] + log1p(exp(-u[large]))
  out[!large] <- log1p(exp(u[!large]))
  out
}

loss_logit <- function(beta, y, X) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  beta <- as.numeric(beta)

  ## beta always contains an intercept slot in position 1
  if (ncol(X) == length(beta) - 1) {
    X <- cbind(1, X)
  }

  if (nrow(X) != length(y)) {
    stop("X and y dimension mismatch in loss_logit")
  }
  if (ncol(X) != length(beta)) {
    stop(sprintf(
      "loss_logit: ncol(X) = %d, length(beta) = %d",
      ncol(X), length(beta)
    ))
  }

  eta <- as.vector(X %*% beta)
  psi_vals <- psi_logit(eta)
  loss <- -2 * sum(y * eta - psi_vals) / length(y)

  if (!is.finite(loss)) {
    return(Inf)
  }

  loss
}

rho_logit <- function(u) {
  rep(1, length(u))
}

loss_logit_compare <- function(beta, X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  beta <- as.numeric(beta)

  if (ncol(X) == length(beta) - 1) {
    X <- cbind(1, X)
  }

  eta <- as.vector(X %*% beta)
  loss <- sum(psi_logit(eta) - y * eta)
  loss / length(y)
}

.drop_na_xy <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)

  keep <- stats::complete.cases(X) & is.finite(y)
  X <- X[keep, , drop = FALSE]
  y <- y[keep]

  if (nrow(X) == 0) {
    stop("All observations were removed after NA filtering")
  }

  list(X = X, y = y)
}

.safe_scale <- function(X) {
  X <- as.matrix(X)
  sds <- apply(X, 2, sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds == 0] <- 1

  mus <- colMeans(X, na.rm = TRUE)
  Xs <- sweep(X, 2, mus, "-")
  Xs <- sweep(Xs, 2, sds, "/")
  Xs[!is.finite(Xs)] <- 0
  Xs
}

.make_lambda_seq <- function(X, y, nlam = 40) {
  X <- as.matrix(X)
  y <- as.numeric(y)

  if (nrow(X) != length(y)) {
    stop("X and y dimension mismatch in make_lambda_seq")
  }

  tmp <- .drop_na_xy(X, y)
  X <- tmp$X
  y <- tmp$y

  Xs <- .safe_scale(X)
  r <- y - mean(y)
  r[!is.finite(r)] <- 0

  grad <- drop(crossprod(Xs, r)) / nrow(Xs)
  grad[!is.finite(grad)] <- 0

  lambda_max <- max(abs(grad), na.rm = TRUE)

  if (!is.finite(lambda_max) || lambda_max <= 0) {
    lambda_max <- 1e-3
  }

  lambda_min_ratio <- 1e-4

  if (nlam <= 1) {
    return(lambda_max)
  }

  lambda_max * (lambda_min_ratio)^((seq_len(nlam) - 1) / (nlam - 1))
}

.safe_mean <- function(x) {
  if (all(is.na(x))) return(Inf)
  mean(x, na.rm = TRUE)
}

transfer_binary_fix <- function(lambda = 1,
                                X, y,
                                t0 = NULL,
                                beta0 = rep(0, ncol(X)),
                                tol = 1e-6,
                                maxit = 400,
                                c = rep(1, ncol(X)),
                                r = 10,
                                use_linesearch = TRUE) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  beta0 <- as.numeric(beta0)
  c <- as.numeric(c)

  p <- ncol(X)

  if (length(y) != nrow(X) || length(beta0) != p || length(c) != p) {
    stop("Dimension mismatch in transfer_binary_fix")
  }

  if (length(lambda) != 1 || is.na(lambda) || !is.finite(lambda) || lambda < 0) {
    stop("lambda must be a finite scalar >= 0")
  }

  beta_prev <- beta0
  y_prev <- beta_prev

  if (is.null(t0)) {
    smax <- base::svd(X, nu = 0, nv = 0)$d[1]
    L <- (smax^2) / (4 * nrow(X))
    t_k <- 1 / L
  } else {
    t_k <- t0
  }

  if (!is.finite(t_k) || t_k <= 0) {
    t_k <- 1
  }

  for (k in 1:maxit) {
    if (use_linesearch) {
      t_try <- t_k
      gy <- grad_g_cpp(y_prev, X, y)
      g_y <- g_cpp(y_prev, X, y)
      repeat {
        G <- G_t_cpp(y_prev, t_try, X, y, lambda, c)
        left <- g_cpp(y_prev - t_try * G, X, y)
        rhs  <- g_y - t_try * sum(gy * G) + (t_try / 2) * sum(G * G)
        if (left <= rhs || t_try < 1e-16) break
        t_try <- 0.5 * t_try
      }
      t_k <- t_try
    }

    gy <- grad_g_cpp(y_prev, X, y)
    v  <- y_prev - t_k * gy
    u  <- soft_thresh_cpp(v, t_k * lambda)
    beta_k <- proj_ker_ct_cpp(u, c)

    y_k <- beta_k + ((k - 1) / (k + r - 1)) * (beta_k - beta_prev)

    if (max(abs(beta_k - beta_prev)) < tol) {
      return(list(beta = beta_k, y = y_k, t = t_k, iter = k, converged = TRUE))
    }

    beta_prev <- beta_k
    y_prev <- y_k
  }

  list(beta = beta_prev, y = y_prev, t = t_k, iter = maxit, converged = FALSE)
}

debias_binary_fix <- function(lambda = 1,
                              X, y,
                              t0 = NULL,
                              beta0 = rep(0, ncol(X)),
                              betaA = rep(0, ncol(X)),
                              tol = 1e-6,
                              maxit = 400,
                              c = rep(1, ncol(X)),
                              r = 10,
                              use_linesearch = TRUE) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  beta0 <- as.numeric(beta0)
  betaA <- as.numeric(betaA)
  c <- as.numeric(c)

  p <- ncol(X)

  if (length(y) != nrow(X) || length(beta0) != p || length(betaA) != p || length(c) != p) {
    stop("Dimension mismatch in debias_binary_fix")
  }

  if (length(lambda) != 1 || is.na(lambda) || !is.finite(lambda) || lambda < 0) {
    stop("lambda must be a finite scalar >= 0")
  }

  beta_prev <- beta0
  y_prev <- beta_prev

  if (is.null(t0)) {
    smax <- base::svd(X, nu = 0, nv = 0)$d[1]
    L <- (smax^2) / (4 * nrow(X))
    t_k <- 1 / L
  } else {
    t_k <- t0
  }

  if (!is.finite(t_k) || t_k <= 0) {
    t_k <- 1
  }

  for (k in 1:maxit) {
    if (use_linesearch) {
      t_try <- t_k
      gy <- grad_g_cpp(y_prev, X, y)
      g_y <- g_cpp(y_prev, X, y)
      repeat {
        G <- G_t_debias_cpp(y_prev, t_try, X, y, lambda, c, betaA)
        left <- g_cpp(y_prev - t_try * G, X, y)
        rhs  <- g_y - t_try * sum(gy * G) + (t_try / 2) * sum(G * G)
        if (left <= rhs || t_try < 1e-16) break
        t_try <- 0.5 * t_try
      }
      t_k <- t_try
    }

    gy <- grad_g_cpp(y_prev, X, y)
    v  <- y_prev - t_k * gy
    u  <- soft_thresh_debias_cpp(v, betaA, t_k * lambda)
    beta_k <- proj_ker_ct_cpp(u, c)

    y_k <- beta_k + ((k - 1) / (k + r - 1)) * (beta_k - beta_prev)

    if (max(abs(beta_k - beta_prev)) < tol) {
      return(list(beta = beta_k,
                  y = y_k,
                  t = t_k,
                  iter = k,
                  converged = TRUE))
    }

    beta_prev <- beta_k
    y_prev <- y_k
  }

  list(beta = beta_prev,
       y = y_prev,
       t = t_k,
       iter = maxit,
       converged = FALSE)
}

transfer_logistic <- function(source = NULL,
                              target = NULL,
                              lambda_beta = NULL, lambda_delta = NULL,
                              nfold = 3, beta_start = NULL, delta_start = NULL, maxit = 600, tol_transfer = 1e-6,
                              tol_debias = 1e-6, Ncov = 0, nlam = 40, intercept = 1, C = NULL) {

  target$x <- as.matrix(target$x)
  target$y <- as.numeric(target$y)

  tmp <- .drop_na_xy(target$x, target$y)
  target$x <- tmp$X
  target$y <- tmp$y

  for (j in seq_along(source)) {
    source[[j]]$x <- as.matrix(source[[j]]$x)
    source[[j]]$y <- as.numeric(source[[j]]$y)
    tmp <- .drop_na_xy(source[[j]]$x, source[[j]]$y)
    source[[j]]$x <- tmp$X
    source[[j]]$y <- tmp$y
  }

  if (is.null(C)) {
    C <- rep(1, times = ncol(target$x) - Ncov)
  }

  if (!intercept) {
    for (i in 0:length(source)) {
      if (i == 0) {
        X_all <- target$x
        y_all <- target$y
      }
      if (i > 0) {
        X_all <- rbind(X_all, source[[i]]$x)
        y_all <- c(y_all, source[[i]]$y)
      }
    }

    reor <- sample(c(1:nrow(X_all)), size = nrow(X_all), replace = FALSE)
    X_all <- X_all[reor, , drop = FALSE]
    y_all <- y_all[reor]
    X_target <- target$x
    y_target <- target$y
    p <- ncol(X_all)

    if (is.null(beta_start)) {
      beta_start <- rep(0, times = p)
    }
    if (is.null(delta_start)) {
      delta_start <- rep(0, times = p)
    }

    Id <- diag(1, nrow = p, ncol = p)
    c_use <- c(rep(0, times = Ncov), C)
    Pc <- c_use %*% t(c_use) / sum(c_use^2)
    X_all <- X_all %*% (Id - Pc)
    X_target <- X_target %*% (Id - Pc)

    tmp_all <- .drop_na_xy(X_all, y_all)
    X_all <- tmp_all$X
    y_all <- tmp_all$y

    tmp_tar <- .drop_na_xy(X_target, y_target)
    X_target <- tmp_tar$X
    y_target <- tmp_tar$y

    if (is.null(lambda_beta)) {
      lambda_beta <- .make_lambda_seq(X_all, y_all, nlam = nlam)
    }
    if (is.null(lambda_delta)) {
      lambda_delta <- .make_lambda_seq(X_target, y_target, nlam = nlam)
    }
  }

  if (intercept) {
    for (i in 0:length(source)) {
      if (i == 0) {
        X_all <- cbind(1, target$x)
        y_all <- target$y
      }
      if (i > 0) {
        X_all <- rbind(X_all, cbind(1, source[[i]]$x))
        y_all <- c(y_all, source[[i]]$y)
      }
    }

    reor <- sample(c(1:nrow(X_all)), size = nrow(X_all), replace = FALSE)
    X_all <- X_all[reor, , drop = FALSE]
    y_all <- y_all[reor]
    X_target <- cbind(1, target$x)
    y_target <- target$y
    p <- ncol(X_all)

    if (is.null(beta_start)) {
      beta_start <- rep(0, times = p)
    }
    if (is.null(delta_start)) {
      delta_start <- rep(0, times = p)
    }

    Id <- diag(1, nrow = p, ncol = p)
    c_use <- c(rep(0, times = Ncov + 1), C)
    Pc <- c_use %*% t(c_use) / sum(c_use^2)
    X_all <- X_all %*% (Id - Pc)
    X_target <- X_target %*% (Id - Pc)

    tmp_all <- .drop_na_xy(X_all, y_all)
    X_all <- tmp_all$X
    y_all <- tmp_all$y

    tmp_tar <- .drop_na_xy(X_target, y_target)
    X_target <- tmp_tar$X
    y_target <- tmp_tar$y

    if (is.null(lambda_beta)) {
      lambda_beta <- .make_lambda_seq(X_all[, -1, drop = FALSE], y_all, nlam = nlam)
    }
    if (is.null(lambda_delta)) {
      lambda_delta <- .make_lambda_seq(X_target[, -1, drop = FALSE], y_target, nlam = nlam)
    }
  }

  loss_trans <- rep(NA_real_, times = length(lambda_beta))
  loss_debias <- rep(NA_real_, times = length(lambda_delta))
  betaA_total <- matrix(0, nrow = length(lambda_beta), ncol = p)

  for (i1 in 1:length(lambda_beta)) {
    loss <- rep(NA_real_, times = nfold)
    lambda_beta1 <- lambda_beta[i1]
    betaA <- NULL

    for (i2 in 1:nfold) {
      c1 <- 1:nrow(X_all)
      m1 <- round(nrow(X_all) / nfold, digits = 0)
      test <- c1[(1 + (i2 - 1) * m1):min((i2 * m1), nrow(X_all))]

      X_all_train <- X_all[-test, , drop = FALSE]
      y_all_train <- y_all[-test]
      X_all_test <- X_all[test, , drop = FALSE]
      y_all_test <- y_all[test]

      fit_try <- tryCatch(
        transfer_binary_fix(lambda = lambda_beta1,
                            X = X_all_train, y = y_all_train,
                            t0 = NULL,
                            beta0 = beta_start,
                            tol = tol_transfer,
                            maxit = maxit,
                            c = c_use,
                            r = 10,
                            use_linesearch = TRUE),
        error = function(e) NULL
      )

      if (!is.null(fit_try)) {
        betaA <- fit_try$beta
        loss[i2] <- loss_logit(beta = betaA, y = y_all_test, X = X_all_test)
      }
    }

    loss_trans[i1] <- .safe_mean(loss)

    if (!is.null(betaA)) {
      betaA_total[i1, ] <- betaA
    }
  }

  good_idx_beta <- which(is.finite(loss_trans))
  if (length(good_idx_beta) == 0) {
    stop("All transfer lambda values failed")
  }

  min_loss_beta <- min(round(loss_trans[good_idx_beta], digits = 4))
  lambda_beta_use <- max(lambda_beta[round(loss_trans, digits = 4) == min_loss_beta])

  betaA <- transfer_binary_fix(lambda = lambda_beta_use,
                               X = X_all, y = y_all,
                               t0 = NULL,
                               beta0 = beta_start,
                               tol = tol_transfer,
                               maxit = maxit,
                               c = c_use,
                               r = 10,
                               use_linesearch = TRUE)
  betaA <- betaA$beta

  delta_total <- matrix(0, nrow = length(lambda_delta), ncol = p)

  for (i1 in 1:length(lambda_delta)) {
    loss <- rep(NA_real_, times = nfold)
    lambda_delta1 <- lambda_delta[i1]
    delta <- NULL

    for (i2 in 1:nfold) {
      c1 <- 1:nrow(X_target)
      m1 <- round(nrow(X_target) / nfold, digits = 0)
      test <- c1[(1 + (i2 - 1) * m1):min((i2 * m1), nrow(X_target))]

      y_target_train <- y_target[-test]
      X_target_train <- X_target[-test, , drop = FALSE]
      y_target_test <- y_target[test]
      X_target_test <- X_target[test, , drop = FALSE]

      fit_try <- tryCatch(
        debias_binary_fix(lambda = lambda_delta1,
                          X = X_target_train, y = y_target_train,
                          t0 = NULL,
                          beta0 = delta_start,
                          betaA = betaA,
                          tol = tol_debias,
                          maxit = maxit,
                          c = c_use,
                          r = 10,
                          use_linesearch = TRUE),
        error = function(e) NULL
      )

      if (!is.null(fit_try)) {
        delta <- fit_try$beta
        loss[i2] <- loss_logit(beta = delta, y = y_target_test, X = X_target_test)
      }
    }

    loss_debias[i1] <- .safe_mean(loss)

    if (!is.null(delta)) {
      delta_total[i1, ] <- delta
    }
  }

  good_idx_delta <- which(is.finite(loss_debias))
  if (length(good_idx_delta) == 0) {
    stop("All debias lambda values failed")
  }

  min_loss_delta <- min(round(loss_debias[good_idx_delta], digits = 4))
  lambda_delta_use <- max(lambda_delta[round(loss_debias, digits = 4) == min_loss_delta])

  delta <- debias_binary_fix(lambda = lambda_delta_use,
                             X = X_target, y = y_target,
                             t0 = NULL,
                             beta0 = delta_start,
                             betaA = betaA,
                             tol = tol_debias,
                             maxit = maxit,
                             c = c_use,
                             r = 10,
                             use_linesearch = TRUE)
  beta_hat <- delta$beta

  if (!intercept) {
    return(list(beta_hat = c(0, beta_hat),
                betaA_total = betaA_total,
                loss_trans = loss_trans,
                loss_debias = loss_debias,
                betaA = c(0, betaA),
                lambda_beta = lambda_beta_use,
                lambda_delta = lambda_delta_use))
  }
  if (intercept) {
    return(list(beta_hat = beta_hat,
                betaA_total = betaA_total,
                loss_trans = loss_trans,
                loss_debias = loss_debias,
                betaA = betaA,
                lambda_beta = lambda_beta_use,
                lambda_delta = lambda_delta_use))
  }
}

transfer_only_logistic <- function(target = NULL,
                                   lambda_beta = NULL,
                                   nfold = 3, beta_start = NULL, maxit = 300, tol_transfer = 1e-4,
                                   Ncov = 0, nlam = 60, intercept = 1, C = NULL) {

  target$x <- as.matrix(target$x)
  target$y <- as.numeric(target$y)

  tmp <- .drop_na_xy(target$x, target$y)
  target$x <- tmp$X
  target$y <- tmp$y

  if (is.null(C)) {
    C <- rep(1, times = ncol(target$x) - Ncov)
  }

  if (!intercept) {
    X_all <- target$x
    y_all <- target$y
    p <- ncol(X_all)

    if (is.null(beta_start)) {
      beta_start <- rep(0, times = p)
    }

    Id <- diag(1, nrow = p, ncol = p)
    c_use <- c(rep(0, times = Ncov), C)
    Pc <- c_use %*% t(c_use) / sum(c_use^2)
    X_all <- X_all %*% (Id - Pc)

    tmp <- .drop_na_xy(X_all, y_all)
    X_all <- tmp$X
    y_all <- tmp$y

    if (is.null(lambda_beta)) {
      lambda_beta <- .make_lambda_seq(X_all, y_all, nlam = nlam)
    }
  }

  if (intercept) {
    X_all <- cbind(1, target$x)
    y_all <- target$y
    p <- ncol(X_all)

    if (is.null(beta_start)) {
      beta_start <- rep(0, times = p)
    }

    Id <- diag(1, nrow = p, ncol = p)
    c_use <- c(rep(0, times = Ncov + 1), C)
    Pc <- c_use %*% t(c_use) / sum(c_use^2)
    X_all <- X_all %*% (Id - Pc)

    tmp <- .drop_na_xy(X_all, y_all)
    X_all <- tmp$X
    y_all <- tmp$y

    if (is.null(lambda_beta)) {
      lambda_beta <- .make_lambda_seq(X_all[, -1, drop = FALSE], y_all, nlam = nlam)
    }
  }

  loss_trans <- rep(NA_real_, times = length(lambda_beta))
  betaA_total <- matrix(0, nrow = length(lambda_beta), ncol = p)

  for (i1 in 1:length(lambda_beta)) {
    loss <- rep(NA_real_, times = nfold)
    lambda_beta1 <- lambda_beta[i1]
    betaA <- NULL

    for (i2 in 1:nfold) {
      c1 <- 1:nrow(X_all)
      m1 <- round(nrow(X_all) / nfold, digits = 0)
      test <- c1[(1 + (i2 - 1) * m1):min((i2 * m1), nrow(X_all))]

      X_all_train <- X_all[-test, , drop = FALSE]
      y_all_train <- y_all[-test]
      X_all_test <- X_all[test, , drop = FALSE]
      y_all_test <- y_all[test]

      fit_try <- tryCatch(
        transfer_binary_fix(lambda = lambda_beta1,
                            X = X_all_train, y = y_all_train,
                            t0 = NULL,
                            beta0 = beta_start,
                            tol = tol_transfer,
                            maxit = maxit,
                            c = c_use,
                            r = 10,
                            use_linesearch = TRUE),
        error = function(e) NULL
      )

      if (!is.null(fit_try)) {
        betaA <- fit_try$beta
        loss[i2] <- loss_logit(beta = betaA, y = y_all_test, X = X_all_test)
      }
    }

    loss_trans[i1] <- .safe_mean(loss)

    if (!is.null(betaA)) {
      betaA_total[i1, ] <- betaA
    }
  }

  good_idx <- which(is.finite(loss_trans))
  if (length(good_idx) == 0) {
    stop("All lambda values failed in transfer_only_logistic")
  }

  min_loss <- min(round(loss_trans[good_idx], digits = 4))
  lambda_beta_use <- max(lambda_beta[round(loss_trans, digits = 4) == min_loss])

  betaA <- transfer_binary_fix(lambda = lambda_beta_use,
                               X = X_all, y = y_all,
                               t0 = NULL,
                               beta0 = beta_start,
                               tol = tol_transfer,
                               maxit = maxit,
                               c = c_use,
                               r = 10,
                               use_linesearch = TRUE)

  beta_hat <- betaA$beta

  if (!intercept) {
    return(list(beta_hat = c(0, beta_hat)))
  }
  if (intercept) {
    return(list(beta_hat = beta_hat))
  }
}

#' Transfer learning for compositional logistic regression under a zero-sum constraint
#'
#' The model supports non-compositional covariates
#' plus compositional covariates whose effects satisfy a zero-sum constraint.
#'
#' @param target_data Target dataset as a list with components `x` and `y`.
#' @param source_data List of source datasets; each element must be a list with
#'   components `x` and `y`.
#' @param nlam Number of candidate lambda values for transfer and debias steps.
#' @param Ncov Number of non-compositional covariates.
#' @param nfold Number of folds used in cross-validation.
#' @param source_id Either `"auto"` or `"all"`.
#' @param intercept Logical; whether to fit an intercept.
#' @param C Optional linear-constraint vector. When `NULL`, the function uses a
#'   vector with zeros for non-compositional covariates and ones for
#'   compositional covariates.
#' @param lambda_list List with components `lambda_transfer` and
#'   `lambda_debias`. Each may be `NULL` to trigger automatic generation.
#' @param maxit Maximum number of optimization iterations.
#' @param tol_transfer Convergence tolerance for the transfer step.
#' @param tol_debias Convergence tolerance for the debias step.
#' @param beta_start Optional starting vector for the transfer step.
#' @param delta_start Optional starting vector for the debias step.
#' @param C0 Optional threshold for selecting transferable sources. When `NULL`,
#'   it is computed from the data.
#'
#' @return If `source_id = "all"`, a list with `beta_hat` and
#'   `transferrable_id`. If `source_id = "auto"`, returns those objects plus
#'   `diff_source`.
#' @export
CatlGLM_binomial <- function(target_data = NULL, source_data = NULL, nlam = 100, Ncov = 0,
                             nfold = 3, source_id = "auto", intercept = TRUE, C = NULL,
                             lambda_list = list(lambda_transfer = NULL, lambda_debias = NULL), maxit = 600,
                             tol_transfer = 1e-6, tol_debias = 1e-6, beta_start = NULL, delta_start = NULL, C0 = NULL) {
  .check_dataset(target_data, "target_data")
  .check_source_list(source_data, "source_data")
  .check_source_id(source_id)

  if (source_id == "auto") {
    LOSS1 <- matrix(NA_real_, nrow = nfold, ncol = length(source_data))

    LOSS1_lasso <- rep(NA_real_, nfold)
    beta_fold_transfer <- vector("list", nfold)
    for (i in 1:nfold) {
      beta_fold_transfer[[i]] <- vector("list", length(source_data))
    }
    score_cv <- matrix(NA_real_, nrow = nfold, ncol = length(source_data))
    beta_constrained_lasso <- vector("list", nfold)

    est_target <- transfer_only_logistic(
      target = target_data,
      lambda_beta = NULL,
      nfold = nfold,
      beta_start = NULL,
      maxit = 400,
      tol_transfer = 1e-5,
      Ncov = Ncov,
      nlam = 50,
      intercept = intercept,
      C = C
    )

    INFO <- info_logistic_n(y = target_data$y, X = target_data$x, beta = est_target$beta_hat)

    for (i in 1:nfold) {
      fold_size <- ceiling(length(target_data$y) / nfold)
      test <- (1 + (i - 1) * fold_size):min(i * fold_size, length(target_data$y))
      train_x <- target_data$x[-test, , drop = FALSE]
      train_y <- target_data$y[-test]
      test_x <- target_data$x[test, , drop = FALSE]
      test_y <- target_data$y[test]

      for (j in 1:length(source_data)) {
        est_j <- tryCatch(
          transfer_only_logistic(
            target = list(x = rbind(train_x, source_data[[j]]$x), y = c(train_y, source_data[[j]]$y)),
            lambda_beta = NULL,
            nfold = nfold,
            beta_start = NULL,
            maxit = 400,
            tol_transfer = 1e-5,
            Ncov = Ncov,
            nlam = 50,
            intercept = intercept,
            C = C
          ),
          error = function(e) NULL
        )

        if (!is.null(est_j)) {
          beta_fold_transfer[[i]][[j]] <- est_j$beta_hat
          LOSS1[i, j] <- loss_logit(beta = est_j$beta_hat, X = test_x, y = test_y)
        }
      }

      lambda_lasso <- tryCatch(
        transfer_only_logistic(
          target = list(x = train_x, y = train_y),
          lambda_beta = NULL,
          nfold = nfold,
          beta_start = NULL,
          maxit = 400,
          tol_transfer = 1e-5,
          Ncov = Ncov,
          nlam = 50,
          intercept = intercept,
          C = C
        ),
        error = function(e) NULL
      )

      if (!is.null(lambda_lasso)) {
        beta_lasso <- lambda_lasso$beta_hat
        LOSS1_lasso[i] <- loss_logit(y = test_y, X = test_x, beta = beta_lasso)
        beta_constrained_lasso[[i]] <- beta_lasso
      }

      print(i)
    }

    for (i in 1:nrow(score_cv)) {
      for (j in 1:ncol(score_cv)) {
        if (!is.null(beta_fold_transfer[[i]][[j]]) &&
            !is.null(beta_constrained_lasso[[i]]) &&
            is.finite(LOSS1_lasso[i]) &&
            LOSS1_lasso[i] > 0) {
          d <- beta_fold_transfer[[i]][[j]] - beta_constrained_lasso[[i]]
          score_cv[i, j] <- as.numeric(t(d) %*% INFO %*% d / LOSS1_lasso[i])
        }
      }
    }

    LOSS_lasso <- mean(LOSS1_lasso, na.rm = TRUE)
    LOSS <- colMeans(LOSS1, na.rm = TRUE)
    score_all <- colMeans(score_cv, na.rm = TRUE)

    if (is.null(C0)) {
      if (all(is.na(score_all))) {
        C0 <- 0
      } else {
        C0 <- min(score_all, na.rm = TRUE)
      }
    }

    diff_source <- LOSS - (1 + C0) * LOSS_lasso
    good_id <- c(1:length(source_data))[is.finite(diff_source) & diff_source <= 0]
    good_source <- source_data[good_id]

    if (length(good_source) == 0) {
      cat("no source are good\n", append = TRUE)
      stop("Error: no valid IDs found length of 'good_id' is 0.")
      utils::flush.console()
    }

    final <- transfer_logistic(
      source = good_source,
      target = target_data,
      lambda_beta = lambda_list$lambda_transfer,
      lambda_delta = lambda_list$lambda_debias,
      nfold = nfold,
      beta_start = beta_start,
      delta_start = delta_start,
      maxit = maxit,
      tol_transfer = tol_transfer,
      tol_debias = tol_debias,
      Ncov = Ncov,
      nlam = nlam,
      C = C,
      intercept = intercept
    )

    return(list(beta_hat = final$beta_hat,
                transferrable_id = good_id,
                diff_source = diff_source))
  }

  if (source_id == "all") {
    final <- transfer_logistic(
      source = source_data,
      target = target_data,
      lambda_beta = lambda_list$lambda_transfer,
      lambda_delta = lambda_list$lambda_debias,
      nfold = nfold,
      beta_start = beta_start,
      delta_start = delta_start,
      maxit = maxit,
      tol_transfer = tol_transfer,
      tol_debias = tol_debias,
      Ncov = Ncov,
      nlam = nlam,
      C = C,
      intercept = intercept
    )
    return(list(beta_hat = final$beta_hat,
                transferrable_id = c(1:length(source_data))))
  }
}
