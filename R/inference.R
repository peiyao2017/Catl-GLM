#' Confidence intervals for compositional transfer-learning estimates
#'
#' Computes de-biased point estimates and Wald confidence intervals for the
#' compositional transfer-learning estimators from [CatlGLM_linear()] and
#' [CatlGLM_binomial()].
#'
#' @param target Target dataset as a list with components `x` and `y`.
#' @param source List of source datasets; each element must be a list with
#'   components `x` and `y`.
#' @param family Model family: `"gaussian"` or `"binomial"`.
#' @param beta.hat Coefficient vector returned by [CatlGLM_linear()] or
#'   [CatlGLM_binomial()].
#' @param cores Number of cores used for parallel computation.
#' @param level Confidence level.
#' @param intercept Logical; must match the fitted model.
#' @param ncov Number of non-compositional covariates.
#'
#' @return A list containing `b.hat`, `beta.hat`, `CI`, and `var.est`.
#' @export
CatlGLM_inf=function (target=NULL, source = NULL,
                              family = c("gaussian", "binomial"), beta.hat = NULL,
                              nodewise.transfer.source.id = "all", cores = 1,
                              level = 0.95, intercept = TRUE,ncov=0) {
  .check_dataset(target, "target")
  .check_source_list(source, "source")
  .check_beta_hat(beta.hat)
  Pc=(c(rep(0,times=ncov),rep(1,times=ncol(target$x)-ncov))/sqrt(sum(rep(1,times=ncol(target$x)-ncov))))%*%(t( c(rep(0,times=ncov),rep(1,times=ncol(target$x)-ncov))/sqrt(sum(rep(1,times=ncol(target$x)-ncov)))))

  target$x=target$x%*%(diag(1,nrow=nrow(Pc),ncol=ncol(Pc))-Pc)

  for(i in 1:length(source)){
    source[[i]]$x=source[[i]]$x%*%(diag(1,nrow=nrow(Pc),ncol=ncol(Pc))-Pc)
  }
  family <- match.arg(family)
  options(warn = 1)
  if (cores <= 1) {
    warning("Only a single core is used. The calculation can be slow, especially when the dimension is large. Multi-cores are suggested.")
  }
  options(warn = 0)
  registerDoParallel(cores)
  r.level <- level + (1 - level)/2
  j <- 0
  if (!is.null(colnames(beta.hat))) {
    beta.hat.names <- colnames(beta.hat)
  }else if (!is.null(names(beta.hat))) {
    beta.hat.names <- names(beta.hat)
  }else {
    beta.hat.names <- NULL
  }


  beta.hat <- as.vector(beta.hat)
  if (!is.null(source) && (is.string(nodewise.transfer.source.id) &&
                           nodewise.transfer.source.id == "all")) {
    nodewise.transfer.source.id <- 1:length(source)
  }else if (is.null(source) || is.null(nodewise.transfer.source.id)) {
    nodewise.transfer.source.id <- NULL
  }

  D <- list(target = target, source = source)
  for (k in 1:(length(nodewise.transfer.source.id) + 1)) {
    if (k == 1) {
      D$target$x <- as.matrix(D$target$x)
    }
    else {
      D$source[[nodewise.transfer.source.id[k - 1]]]$x <- as.matrix(D$source[[nodewise.transfer.source.id[k -
                                                                                                            1]]]$x)
    }
  }
  if (family == "gaussian") {
    disper=sqrt(sum((target$y-cbind(1,target$x)%*%beta.hat)^2)/nrow(target$x))
    D.centeralized <- D
    if (intercept) {
      for (k in 1:(length(nodewise.transfer.source.id) +
                   1)) {
        if (k > 1) {
          D.centeralized$source[[nodewise.transfer.source.id[k -
                                                               1]]]$x <- cbind(1, D.centeralized$source[[nodewise.transfer.source.id[k -
                                                                                                                                       1]]]$x)
        }
        else {
          D.centeralized$target$x <- cbind(1, D.centeralized$target$x)
        }
      }
    }
    p <- ncol(D.centeralized$target$x)
    X.comb <- foreach(k = unique(c(0, nodewise.transfer.source.id)),
                      .combine = "rbind") %do% {
                        if (k > 0) {
                          D.centeralized$source[[k]]$x
                        }
                        else {
                          D.centeralized$target$x
                        }
                      }
    Sigma.hat <- t(X.comb) %*% X.comb/nrow(X.comb)
    L <- foreach(j = 1:p, .combine = "rbind") %dopar% {
            D1 <- D.centeralized
      for (k in 1:(length(nodewise.transfer.source.id) +
                   1)) {
        if (k > 1) {
          D1$source[[nodewise.transfer.source.id[k -
                                                   1]]]$y <- D1$source[[nodewise.transfer.source.id[k -
                                                                                                      1]]]$x[, j]
          D1$source[[nodewise.transfer.source.id[k -
                                                   1]]]$x <- D1$source[[nodewise.transfer.source.id[k -
                                                                                                      1]]]$x[, -j]
        }
        else {
          D1$target$y <- D1$target$x[, j]
          D1$target$x <- D1$target$x[, -j]
        }
      }
      node.lasso <- glmtrans::glmtrans(target = D1$target, source = D1$source,
                             family = "gaussian", alg = "ori", transfer.source.id = nodewise.transfer.source.id,
                             intercept = FALSE, detection.info = FALSE )
      gamma <- node.lasso$beta[-1]
      tau2 <- Sigma.hat[j, j] - Sigma.hat[j, -j, drop = F] %*%
        gamma
      theta <- rep(1, p)
      theta[-j] <- -gamma
      c(theta, tau2)
    }
    Theta.hat <- solve(diag(L[, p + 1])) %*% L[1:p, 1:p]
    if(intercept){
      Pc_int=(c(rep(0,times=ncov+1),rep(1,times=ncol(target$x)-ncov))/sqrt(sum(rep(1,times=ncol(target$x)-ncov))))%*%(t( c(rep(0,times=ncov+1),rep(1,times=ncol(target$x)-ncov))/sqrt(sum(rep(1,times=ncol(target$x)-ncov)))))
      Theta.tilde=(diag(1,nrow=nrow(Pc_int),ncol=ncol(Pc_int))-Pc_int)%*%Theta.hat%*%(diag(1,nrow=nrow(Pc_int),ncol=ncol(Pc_int))-Pc_int)
    }
    if(!intercept){

      Theta.tilde=(diag(1,nrow=nrow(Pc),ncol=ncol(Pc))-Pc)%*%Theta.hat%*%(diag(1,nrow=nrow(Pc),ncol=ncol(Pc))-Pc)
    }
    Z <- D$target$y - D$target$x %*% beta.hat[-1] - beta.hat[1]
    if (intercept) {
      b.hat <- as.matrix(beta.hat) + Theta.tilde %*% t(cbind(1,
                                                           D$target$x)) %*% Z/length(D$target$y)
    }
    else {
      b.hat <- beta.hat[-1] + Theta.tilde %*% t(D$target$x) %*%
        Z/length(D$target$y)
    }

    var.est <- diag(Theta.tilde %*% Sigma.hat %*% t(Theta.tilde))
    CI <- data.frame(b.hat = b.hat, lb = b.hat - disper*qnorm(r.level) *
                       sqrt(var.est/length(D$target$y)), ub = b.hat + disper*qnorm(r.level) *
                       sqrt(var.est/length(D$target$y)))
  }
  else if (family == "binomial") {
    Dw <- D
    if (intercept) {
      for (k in 1:(length(nodewise.transfer.source.id) +
                   1)) {
        if (k > 1) {
          uk <- Dw$source[[nodewise.transfer.source.id[k -
                                                         1]]]$x %*% beta.hat[-1] + beta.hat[1]
          wk <- as.vector(sqrt(exp(-1*uk)/((1 + exp(-uk))^2)))
          Dw$source[[nodewise.transfer.source.id[k -
                                                   1]]]$x <- cbind(wk, diag(wk, nrow = length(wk)) %*%
                                                                     Dw$source[[nodewise.transfer.source.id[k -
                                                                                                              1]]]$x)
        }
        else {
          uk <- Dw$target$x %*% beta.hat[-1] + beta.hat[1]
          wk <- as.vector(sqrt(exp(-1*uk)/((1 + exp(-uk))^2)))
          Dw$target$x <- cbind(wk, diag(wk, nrow = length(wk)) %*%
                                 Dw$target$x)
        }
      }
    }else {
      for (k in 1:(length(nodewise.transfer.source.id) +
                   1)) {
        if (k > 1) {
          uk <- Dw$source[[nodewise.transfer.source.id[k -
                                                         1]]]$x %*% beta.hat[-1]
          wk <- as.vector(sqrt(exp(-1*uk)/((1 + exp(-uk))^2)))
          Dw$source[[nodewise.transfer.source.id[k -
                                                   1]]]$x <- diag(wk) %*% Dw$source[[nodewise.transfer.source.id[k -
                                                                                                                   1]]]$x
        }
        else {
          uk <- Dw$target$x %*% beta.hat[-1]
          wk <- as.vector(sqrt(exp(-1*uk)/((1 + exp(-uk))^2)))
          Dw$target$x <- diag(wk) %*% Dw$target$x
        }
      }
    }

    Xw <- foreach(k = unique(c(0, nodewise.transfer.source.id)),
                  .combine = "rbind") %do% {
                    if (k > 0) {
                      Dw$source[[k]]$x
                    }
                    else {
                      Dw$target$x
                    }
                  }
    p <- ncol(Dw$target$x)
    Sigma.hat <- (t(Xw) %*% Xw)/nrow(Xw)
    L <- foreach(j = 1:p, .combine = "rbind") %dopar% {
            D1 <- Dw
      for (k in 1:(length(nodewise.transfer.source.id) +
                   1)) {
        if (k > 1) {
          D1$source[[nodewise.transfer.source.id[k -
                                                   1]]]$y <- D1$source[[nodewise.transfer.source.id[k -
                                                                                                      1]]]$x[, j]
          D1$source[[nodewise.transfer.source.id[k -
                                                   1]]]$x <- D1$source[[nodewise.transfer.source.id[k -
                                                                                                      1]]]$x[, -j]
        }
        else {
          D1$target$y <- D1$target$x[, j]
          if (all(abs(D1$target$y) <= 1e-20)) {
            D1$target$y <- rep(1e-20, length(D1$target$y))
          }
          D1$target$x <- D1$target$x[, -j]
        }
      }
      node.lasso <- try(glmtrans::glmtrans(target = D1$target, source = D1$source,
                                 family = "gaussian", alg = "ori", transfer.source.id = nodewise.transfer.source.id,
                                 intercept = FALSE ))
      if (inherits(node.lasso, "try-error")) {
        stop(paste("errors happened in feature", j))
      }
      gamma <- node.lasso$beta[-1]
      tau2 <- Sigma.hat[j, j] - Sigma.hat[j, -j, drop = F] %*%
        gamma
      theta <- rep(1, p)
      theta[-j] <- -gamma
      c(theta, tau2)
    }
    tau2 <- L[, p + 1]
    tau2.inv <- 1/tau2
    tau2.inv[abs(tau2) <= 1e-20] <- 0
    Theta.hat <- diag(tau2.inv) %*% L[1:p, 1:p]
    if(intercept){
      Pc_int=(c(rep(0,times=ncov+1),rep(1,times=ncol(target$x)-ncov))/sqrt(sum(rep(1,times=ncol(target$x)-ncov))))%*%(t( c(rep(0,times=ncov+1),rep(1,times=ncol(target$x)-ncov))/sqrt(sum(rep(1,times=ncol(target$x)-ncov)))))
      Theta.tilde=(diag(1,nrow=nrow(Pc_int),ncol=ncol(Pc_int))-Pc_int)%*%Theta.hat%*%(diag(1,nrow=nrow(Pc_int),ncol=ncol(Pc_int))-Pc_int)
    }
    if(!intercept){
      Theta.tilde=(diag(1,nrow=nrow(Pc),ncol=ncol(Pc))-Pc)%*%Theta.hat%*%(diag(1,nrow=nrow(Pc),ncol=ncol(Pc))-Pc)
    }
    u.target <- D$target$x %*% beta.hat[-1] + beta.hat[1]
    Z <- D$target$y - 1/(1 + exp(-u.target))
    if (intercept) {
      b.hat <- as.matrix(beta.hat) + Theta.tilde %*% t(cbind(1,
                                                           D$target$x)) %*% Z/length(D$target$y)
    }else {
      b.hat <- beta.hat[-1] + Theta.tilde %*% t(D$target$x) %*%
        Z/length(D$target$y)
    }


    var.est <- diag(Theta.tilde %*% Sigma.hat %*% t(Theta.tilde))
    CI <- data.frame(b.hat = b.hat, lb = b.hat - qnorm(r.level) *
                       sqrt(var.est/length(D$target$y)), ub = b.hat + qnorm(r.level) *
                       sqrt(var.est/length(D$target$y)))
  }
  else if (family == "poisson") {
    Dw <- D
    if (intercept) {
      for (k in 1:(length(nodewise.transfer.source.id) +
                   1)) {
        if (k > 1) {
          uk <- Dw$source[[nodewise.transfer.source.id[k -
                                                         1]]]$x %*% beta.hat[-1] + beta.hat[1]
          wk <- as.vector(exp(uk/2))
          Dw$source[[nodewise.transfer.source.id[k -
                                                   1]]]$x <- cbind(wk, diag(wk) %*% Dw$source[[nodewise.transfer.source.id[k -
                                                                                                                             1]]]$x)
        }
        else {
          uk <- Dw$target$x %*% beta.hat[-1] + beta.hat[1]
          wk <- as.vector(exp(uk/2))
          Dw$target$x <- cbind(wk, diag(wk) %*% Dw$target$x)
        }
      }
    }
    Xw <- foreach(k = unique(c(0, nodewise.transfer.source.id)),
                  .combine = "rbind") %do% {
                    if (k > 0) {
                      Dw$source[[k]]$x
                    }
                    else {
                      Dw$target$x
                    }
                  }
    p <- ncol(Dw$target$x)
    Sigma.hat <- (t(Xw) %*% Xw)/nrow(Xw)
    L <- foreach(j = 1:p, .combine = "rbind") %dopar% {
            D1 <- Dw
      for (k in 1:(length(nodewise.transfer.source.id) +
                   1)) {
        if (k > 1) {
          D1$source[[nodewise.transfer.source.id[k -
                                                   1]]]$y <- D1$source[[nodewise.transfer.source.id[k -
                                                                                                      1]]]$x[, j]
          D1$source[[nodewise.transfer.source.id[k -
                                                   1]]]$x <- D1$source[[nodewise.transfer.source.id[k -
                                                                                                      1]]]$x[, -j]
        }
        else {
          D1$target$y <- D1$target$x[, j]
          D1$target$x <- D1$target$x[, -j]
        }
      }
      node.lasso <- glmtrans::glmtrans(target = D1$target, source = D1$source,
                             family = "gaussian", alg = "ori", transfer.source.id = nodewise.transfer.source.id,
                             intercept = FALSE )
      gamma <- node.lasso$beta[-1]
      tau2 <- Sigma.hat[j, j] - Sigma.hat[j, -j, drop = F] %*%
        gamma
      theta <- rep(1, p)
      theta[-j] <- -gamma
      c(theta, tau2)
    }
    Theta.hat <- solve(diag(L[, p + 1])) %*% L[1:p, 1:p]
    if(intercept){
      Pc_int=(c(rep(0,times=ncov+1),rep(1,times=ncol(target$x)-ncov))/sqrt(sum(rep(1,times=ncol(target$x)-ncov))))%*%(t( c(rep(0,times=ncov+1),rep(1,times=ncol(target$x)-ncov))/sqrt(sum(rep(1,times=ncol(target$x)-ncov)))))
      Theta.tilde=(diag(1,nrow=nrow(Pc_int),ncol=ncol(Pc_int))-Pc_int)%*%Theta.hat%*%(diag(1,nrow=nrow(Pc_int),ncol=ncol(Pc_int))-Pc_int)
    }
    if(!intercept){

      Theta.tilde=(diag(1,nrow=nrow(Pc),ncol=ncol(Pc))-Pc)%*%Theta.hat%*%(diag(1,nrow=nrow(Pc),ncol=ncol(Pc))-Pc)
    }
    u.target <- D$target$x %*% beta.hat[-1] + beta.hat[1]
    Z <- D$target$y - exp(u.target)
    if (intercept) {
      b.hat <- as.matrix(beta.hat) + Theta.tilde %*% t(cbind(1,
                                                           D$target$x)) %*% Z/length(D$target$y)
    }
    else {
      b.hat <- beta.hat[-1] + Theta.tilde%*% t(D$target$x) %*%
        Z/length(D$target$y)
    }

    var.est <- diag(Theta.tilde %*% Sigma.hat %*% t(Theta.tilde))
    CI <- data.frame(b.hat = b.hat, lb = b.hat - qnorm(r.level) *
                       sqrt(var.est/length(D$target$y)), ub = b.hat + qnorm(r.level) *
                       sqrt(var.est/length(D$target$y)))
  }
  stopImplicitCluster()
  if (!is.null(beta.hat.names)) {
    rownames(CI) <- beta.hat.names
    names(var.est) <- beta.hat.names
    names(b.hat) <- beta.hat.names
    names(beta.hat) <- beta.hat.names
  }
  return(list(b.hat = b.hat, beta.hat = beta.hat, CI = CI,
              var.est = var.est))
}
