
info_linear_n <- function(y, X, beta, sigma2 = NULL,intercept=TRUE) {
  if(intercept){
    X=cbind(1,X)
  }
  X <- as.matrix(X)
  y <- as.numeric(y)
  beta <- as.numeric(beta)

  n <- length(y)

  if (nrow(X) != n)
    stop("X and y dimension mismatch")
  if (ncol(X) != length(beta))
    stop("beta length mismatch")

  mu <- as.vector(X %*% beta)

  # variance estimate (MLE-type; works when p > n)
  if (is.null(sigma2))
    sigma2 <- sum((y - mu)^2) / n

  # normalized information
  I <- crossprod(X) / (n * sigma2)

  return(I)
}

loss_linear=function(beta,X,y){
  a=2*sum((y-X%*%beta)^2)/length(y)
  return(a)
}

## ================================================================
## 1) C++ for transfer_linear_fix (shared helpers + G_t_ls_cpp)
## ================================================================


## ================================================================
## 2) R function: transfer_linear_fix
## ================================================================
transfer_linear_fix <- function(lambda = 1,
                                X, y,
                                t0 = NULL,
                                beta0 = rep(0, ncol(X)),
                                tol = 1e-6,
                                maxit = 400,
                                c = rep(1, ncol(X)),
                                r = 10,
                                use_linesearch = TRUE,
                                monotone_restart = TRUE,
                                return_trace = TRUE) {
  X <- as.matrix(X); y <- as.numeric(y)
  beta0 <- as.numeric(beta0); c <- as.numeric(c)
  n <- nrow(X); p <- ncol(X)
  stopifnot(length(y)==n, length(beta0)==p, length(c)==p)
  if (lambda < 0) stop("lambda must be >= 0")
  
  beta_prev <- beta0
  y_prev    <- beta0
  
  # Lipschitz: L = ||X||_op^2 / n
  t_k <- if (is.null(t0)) {
    smax <- base::svd(X, nu = 0, nv = 0)$d[1]; 1 / ((smax^2) / n)
  } else t0
  
  if (return_trace) {
    obj_hist <- numeric(maxit); step_hist <- numeric(maxit); gap_hist <- numeric(maxit)
  }
  g_y_prev <- g_ls_cpp(y_prev, X, y)
  
  for (k in 1:maxit) {
    if (use_linesearch) {
      t_try <- t_k
      gy <- grad_g_ls_cpp(y_prev, X, y)
      repeat {
        G <- G_t_ls_constrained_cpp(y_prev, t_try, X, y, lambda, c)
        left <- g_ls_cpp(y_prev - t_try * G, X, y)
        rhs  <- g_y_prev - t_try * sum(gy * G) + (t_try/2) * sum(G*G)
        if (left <= rhs || t_try < 1e-16) break
        t_try <- 0.5 * t_try
      }
      t_k <- t_try
    } else {
      gy <- grad_g_ls_cpp(y_prev, X, y)
    }
    
    v      <- y_prev - t_k * gy
    beta_k <- prox_l1_zerosum_rank1_cpp(v, t_k * lambda, c)
    
    if (monotone_restart) {
      g_beta_k <- g_ls_cpp(beta_k, X, y)
      if (g_beta_k > g_y_prev) { y_prev <- beta_prev; g_y_prev <- g_ls_cpp(y_prev, X, y); next }
      g_y_prev <- g_beta_k
    } else {
      g_y_prev <- g_ls_cpp(beta_k, X, y)
    }
    
    y_k <- beta_k + ((k - 1)/(k + r - 1)) * (beta_k - beta_prev)
    
    if (max(abs(beta_k - beta_prev)) < tol) {
      out <- list(beta = beta_k, y = y_k, t = t_k, iter = k, converged = TRUE)
      if (return_trace) {
        obj_hist[k] <- g_y_prev; step_hist[k] <- t_k; gap_hist[k] <- max(abs(beta_k - beta_prev))
        out$trace <- list(obj = obj_hist[1:k], step = step_hist[1:k], max_change = gap_hist[1:k])
      }
      return(out)
    }
    
    if (return_trace) {
      obj_hist[k] <- g_y_prev; step_hist[k] <- t_k; gap_hist[k] <- max(abs(beta_k - beta_prev))
    }
    beta_prev <- beta_k; y_prev <- y_k
  }
  
  out <- list(beta = beta_prev, y = y_prev, t = t_k, iter = maxit, converged = FALSE)
  if (return_trace) out$trace <- list(obj = obj_hist, step = step_hist, max_change = gap_hist)
  out
}


## ================================================================
## 3) C++ for transfer_linear_debias_fix (debias prox + mapping)
## ================================================================


## ================================================================
## 4) R function: transfer_linear_debias_fix
## ================================================================
debias_linear_fix <- function(lambda = 1,
                              X, y,
                              betaA,                 # anchor vector
                              t0 = NULL,
                              beta0 = rep(0, ncol(X)),
                              tol = 1e-6,
                              maxit = 400,
                              c = rep(1, ncol(X)),   # 1 = constrained; 0 = unconstrained
                              r = 10,
                              use_linesearch = TRUE,
                              monotone_restart = TRUE,
                              return_trace = TRUE) {
  X <- as.matrix(X); y <- as.numeric(y)
  beta0 <- as.numeric(beta0); c <- as.numeric(c); betaA <- as.numeric(betaA)
  n <- nrow(X); p <- ncol(X)
  stopifnot(length(y)==n, length(beta0)==p, length(c)==p, length(betaA)==p)
  if (lambda < 0) stop("lambda must be >= 0")
  
  beta_prev <- beta0
  y_prev    <- beta0
  
  # Lipschitz: L = ||X||_op^2 / n
  t_k <- if (is.null(t0)) {
    smax <- base::svd(X, nu = 0, nv = 0)$d[1]; 1 / ((smax^2) / n)
  } else t0
  
  if (return_trace) {
    obj_hist <- numeric(maxit); step_hist <- numeric(maxit); gap_hist <- numeric(maxit)
  }
  g_y_prev <- g_ls_cpp(y_prev, X, y)
  
  for (k in 1:maxit) {
    # ---- backtracking line search (debias mapping) ----
    if (use_linesearch) {
      t_try <- t_k
      gy <- grad_g_ls_cpp(y_prev, X, y)
      repeat {
        G <- G_t_ls_debias_constrained_cpp(y_prev, t_try, X, y, lambda, c, betaA)
        left <- g_ls_cpp(y_prev - t_try * G, X, y)
        rhs  <- g_y_prev - t_try * sum(gy * G) + (t_try / 2) * sum(G * G)
        if (left <= rhs || t_try < 1e-16) break
        t_try <- 0.5 * t_try
      }
      t_k <- t_try
    } else {
      gy <- grad_g_ls_cpp(y_prev, X, y)
    }
    
    # ---- proximal update using exact anchored constrained prox ----
    v       <- y_prev - t_k * gy
    beta_k  <- prox_l1_shift_zerosum_rank1_cpp(v, t_k * lambda, c, betaA)
    
    # ---- optional monotone restart ----
    if (monotone_restart) {
      g_beta_k <- g_ls_cpp(beta_k, X, y)
      if (g_beta_k > g_y_prev) {
        y_prev   <- beta_prev
        g_y_prev <- g_ls_cpp(y_prev, X, y)
        next
      } else {
        g_y_prev <- g_beta_k
      }
    } else {
      g_y_prev <- g_ls_cpp(beta_k, X, y)
    }
    
    # ---- Nesterov acceleration ----
    y_k <- beta_k + ((k - 1) / (k + r - 1)) * (beta_k - beta_prev)
    
    # ---- stopping ----
    if (max(abs(beta_k - beta_prev)) < tol) {
      out <- list(beta = beta_k, y = y_k, t = t_k, iter = k, converged = TRUE)
      if (return_trace) {
        obj_hist[k]  <- g_y_prev
        step_hist[k] <- t_k
        gap_hist[k]  <- max(abs(beta_k - beta_prev))
        out$trace <- list(obj = obj_hist[1:k], step = step_hist[1:k], max_change = gap_hist[1:k])
      }
      return(out)
    }
    
    if (return_trace) {
      obj_hist[k]  <- g_y_prev
      step_hist[k] <- t_k
      gap_hist[k]  <- max(abs(beta_k - beta_prev))
    }
    beta_prev <- beta_k
    y_prev    <- y_k
  }
  
  out <- list(beta = beta_prev, y = y_prev, t = t_k, iter = maxit, converged = FALSE)
  if (return_trace) out$trace <- list(obj = obj_hist, step = step_hist, max_change = gap_hist)
  out
}






transfer_linear=function(source=NULL,
                           target=NULL,
                           lambda_beta=NULL,lambda_delta=NULL,
                           nfold=3,beta_start=NULL,delta_start=NULL,maxit=600,tol_transfer=1e-6,
                           tol_debias=1e-6,Ncov=0,nlam=40,intercept=1,C=NULL){
  
  
  if(is.null(C)){
    C=rep(1,times=ncol(target$x)-Ncov)
  }
  if(!intercept){
    for(i in 0:length(source)){
      if(i==0){
        X_all=target$x
        y_all=target$y
      }
      if(i>0){
        X_all=rbind(X_all,source[[i]]$x)
        y_all=c(y_all,source[[i]]$y)
      }
    }
    reor=sample(c(1:nrow(X_all)),size=nrow(X_all),replace = FALSE)
    X_all=X_all[reor,]
    y_all=y_all[reor]
    X_target=target$x
    y_target=target$y
    p=ncol(X_all)
    Nmicro=p-Ncov
    
    
    if(is.null(beta_start)){
      beta_start=rep(0,times=p)
    }
    if(is.null(delta_start)){
      delta_start=rep(0,times=p)
    }
    Id=diag(1,nrow=p,ncol=p)
    c_use=c(rep(0,times=Ncov),C)
    Pc=c_use%*%t(c_use)/Nmicro
    X_all=X_all%*%(Id-Pc)
    X_target=X_target%*%(Id-Pc)
    
    
    if(is.null(lambda_beta)){
      x_std_transfer=base::scale(X_all,center = TRUE,scale=TRUE)
      p_hat_transfer=mean(y_all)
      r_transfer=y_all - p_hat_transfer
      alpha=1
      grad=drop(crossprod(x_std_transfer, r_transfer)) / nrow(x_std_transfer)
      lambda_max_transfer=max(abs(grad)) / alpha
      lambda_min_ratio_transfer=if (nrow(x_std_transfer) > p) 1e-4 else 1e-4
      lambda_beta=lambda_max_transfer * (lambda_min_ratio_transfer) ^ ((seq_len(nlam) - 1) / (nlam - 1))
    }
    if(is.null(lambda_delta)){
      x_std_debias=base::scale(X_target,center = TRUE,scale=TRUE)
      p_hat_debias=mean(y_target)
      r_debias=y_target - p_hat_debias
      alpha=1
      grad=drop(crossprod(x_std_debias, r_debias)) / nrow(x_std_debias)
      lambda_max_debias=max(abs(grad)) / alpha
      lambda_min_ratio_debias=if (nrow(x_std_debias) > p) 1e-4 else 1e-4
      lambda_delta=lambda_max_debias * (lambda_min_ratio_debias) ^ ((seq_len(nlam) - 1) / (nlam - 1))
    }
  }
  
  if(intercept){
    for(i in 0:length(source)){
      if(i==0){
        X_all=cbind(1,target$x)
        y_all=target$y
      }
      if(i>0){
        X_all=rbind(X_all,cbind(1,source[[i]]$x))
        y_all=c(y_all,source[[i]]$y)
      }
    }
    reor=sample(c(1:nrow(X_all)),size=nrow(X_all),replace = FALSE)
    X_all=X_all[reor,]
    y_all=y_all[reor]
    X_target=cbind(1,target$x)
    y_target=target$y
    p=ncol(X_all)
    Nmicro=p-Ncov-1
    
    
    if(is.null(beta_start)){
      beta_start=rep(0,times=p)
    }
    if(is.null(delta_start)){
      delta_start=rep(0,times=p)
    }
    Id=diag(1,nrow=p,ncol=p)
    c_use=c(rep(0,times=Ncov+1),C)
    Pc=c_use%*%t(c_use)/Nmicro
    X_all=X_all%*%(Id-Pc)
    X_target=X_target%*%(Id-Pc)
    
    if(is.null(lambda_beta)){
      x_std_transfer=base::scale(X_all[,-1],center = TRUE,scale=TRUE)
      p_hat_transfer=mean(y_all)
      r_transfer=y_all - p_hat_transfer
      alpha=1
      grad=drop(crossprod(x_std_transfer, r_transfer)) / nrow(x_std_transfer)
      lambda_max_transfer=max(abs(grad)) / alpha
      lambda_min_ratio_transfer=if (nrow(x_std_transfer) > p) 1e-4 else 1e-4
      lambda_beta=lambda_max_transfer * (lambda_min_ratio_transfer) ^ ((seq_len(nlam) - 1) / (nlam - 1))
    }
    if(is.null(lambda_delta)){
      x_std_debias=base::scale(X_target[,-1],center = TRUE,scale=TRUE)
      p_hat_debias=mean(y_target)
      r_debias=y_target - p_hat_debias
      alpha=1
      grad=drop(crossprod(x_std_debias, r_debias)) / nrow(x_std_debias)
      lambda_max_debias=max(abs(grad)) / alpha
      lambda_min_ratio_debias=if (nrow(x_std_debias) > p) 1e-4 else 1e-4
      lambda_delta=lambda_max_debias * (lambda_min_ratio_debias) ^ ((seq_len(nlam) - 1) / (nlam - 1))
    }
  }
  
  loss_trans=rep(0,times=length(lambda_beta))
  loss_debias=rep(0,times=length(lambda_delta))
  
  betaA_total=matrix(0,nrow=length(lambda_beta),ncol=p)
  
  for(i1 in 1:length(lambda_beta)){
    
    loss=rep(0,times=nfold)
    lambda_beta1=lambda_beta[i1]
    for(i2 in 1:nfold){
      
      c1=1:nrow(X_all)
      m1=round(nrow(X_all)/nfold,digits = 0)
      test=c1[(1+(i2-1)*m1):min((i2*m1),nrow(X_all))]
      X_all_train=X_all[-test,]
      y_all_train=y_all[-test]
      X_all_test=X_all[test,]
      y_all_test=y_all[test] 
      
      betaA=transfer_linear_fix(lambda = lambda_beta1,
                                X=X_all_train, y=y_all_train,
                                t0 = NULL,
                                beta0 = beta_start,
                                tol = tol_transfer,
                                maxit = maxit,
                                c = c_use,
                                r = 10,
                                use_linesearch = TRUE)
      betaA=betaA$beta
      loss[i2]=loss_linear(beta=betaA,y=y_all_test,X=X_all_test) 
    }
    loss_trans[i1]=mean(loss)
    
    betaA_total[i1,]=betaA
    
  }
  lambda_beta_use=max(lambda_beta[round(loss_trans,digits = 4)==min(round(loss_trans,digits = 4))])
  lambda_beta_index=max(c(1:length(lambda_beta))[round(loss_trans,digits = 4)==min(round(loss_trans,digits = 4))])
  betaA=transfer_linear_fix(lambda = lambda_beta_use,
                            X=X_all, y=y_all,
                            t0 = NULL,
                            beta0 = beta_start,
                            tol = tol_transfer,
                            maxit = maxit,
                            c = c_use,
                            r = 10,
                            use_linesearch = TRUE)
  betaA=betaA$beta
  #betaA=betaA_total[lambda_beta_index,]
  delta_total=matrix(0,nrow=length(lambda_delta),ncol=p)
  
  for(i1 in 1:length(lambda_delta)){
    
    loss=rep(0,times=nfold)
    lambda_delta1=lambda_delta[i1]
    for(i2 in 1:nfold){
      
      c1=1:nrow(X_target)
      m1=round(nrow(X_target)/nfold,digits = 0)
      test=c1[(1+(i2-1)*m1):min((i2*m1),nrow(X_target))]
      y_target_train=y_target[-test]
      X_target_train=X_target[-test,]
      y_target_test=y_target[test]
      X_target_test=X_target[test,]
      
      
      
      delta=debias_linear_fix(lambda = lambda_delta1,
                              X=X_target_train, y=y_target_train,
                              t0 = NULL,
                              beta0 = delta_start,
                              betaA = betaA,
                              tol = tol_debias,
                              maxit = maxit,
                              c = c_use,
                              r = 10,
                              use_linesearch = TRUE)
      delta=delta$beta
      loss[i2]=loss_linear(beta=delta,y=y_target_test,X=X_target_test) 
      
    }
    loss_debias[i1]=mean(loss)
    delta_total[i1,]=delta
    
    
    
  }
  lambda_delta_use=max(lambda_delta[round(loss_debias,digits = 4)==min(round(loss_debias,digits = 4))])
  lambda_delta_index=max(c(1:length(lambda_delta))[round(loss_debias,digits = 4)==min(round(loss_debias,digits = 4))])
  delta=debias_linear_fix(lambda = lambda_delta_use,
                          X=X_target, y=y_target,
                          t0 = NULL,
                          beta0 = delta_start,
                          betaA = betaA,
                          tol = tol_debias,
                          maxit = maxit,
                          c = c_use,
                          r = 10,
                          use_linesearch = TRUE)
  beta_hat=delta$beta
  #beta_hat=delta_total[lambda_delta_index,]
  
  if(!intercept){
    return(list(beta_hat=c(0,beta_hat),betaA_total=betaA_total,loss_trans=loss_trans,loss_debias=loss_debias,betaA=c(0,betaA),lambda_beta=lambda_beta_use,lambda_delta=lambda_delta_use))
  }
  if(intercept){
    return(list(beta_hat=beta_hat,betaA_total=betaA_total,loss_trans=loss_trans,loss_debias=loss_debias,betaA=betaA,lambda_beta=lambda_beta_use,lambda_delta=lambda_delta_use))
  }
  
  
}




transfer_only_linear=function( 
    target=NULL,
    lambda_beta=NULL, 
    nfold=3,beta_start=NULL,maxit=300,tol_transfer=1e-4,
    Ncov=0,nlam=60,intercept=1,C=NULL){
  
  
  if(is.null(C)){
    C=rep(1,times=ncol(target$x)-Ncov)
  }
  if(!intercept){
    
    X_all=target$x
    y_all=target$y
    
    
    
    p=ncol(X_all)
    Nmicro=p-Ncov
    
    
    if(is.null(beta_start)){
      beta_start=rep(0,times=p)
    }
    
    Id=diag(1,nrow=p,ncol=p)
    c_use=c(rep(0,times=Ncov),C)
    Pc=c_use%*%t(c_use)/Nmicro
    X_all=X_all%*%(Id-Pc)
    
    
    if(is.null(lambda_beta)){
      x_std_transfer=base::scale(X_all,center = TRUE,scale=TRUE)
      p_hat_transfer=mean(y_all)
      r_transfer=y_all - p_hat_transfer
      alpha=1
      grad=drop(crossprod(x_std_transfer, r_transfer)) / nrow(x_std_transfer)
      lambda_max_transfer=max(abs(grad)) / alpha
      lambda_min_ratio_transfer=if (nrow(x_std_transfer) > p) 1e-4 else 1e-4
      lambda_beta=lambda_max_transfer * (lambda_min_ratio_transfer) ^ ((seq_len(nlam) - 1) / (nlam - 1))
    }
    
  }
  
  if(intercept){
    
    X_all=cbind(1,target$x)
    y_all=target$y
    
    
    p=ncol(X_all)
    Nmicro=p-Ncov-1
    
    
    if(is.null(beta_start)){
      beta_start=rep(0,times=p)
    }
    
    Id=diag(1,nrow=p,ncol=p)
    c_use=c(rep(0,times=Ncov+1),C)
    Pc=c_use%*%t(c_use)/Nmicro
    X_all=X_all%*%(Id-Pc)
    
    
    if(is.null(lambda_beta)){
      x_std_transfer=base::scale(X_all[,-1],center = TRUE,scale=TRUE)
      p_hat_transfer=mean(y_all)
      r_transfer=y_all - p_hat_transfer
      alpha=1
      grad=drop(crossprod(x_std_transfer, r_transfer)) / nrow(x_std_transfer)
      lambda_max_transfer=max(abs(grad)) / alpha
      lambda_min_ratio_transfer=if (nrow(x_std_transfer) > p) 1e-4 else 1e-4
      lambda_beta=lambda_max_transfer * (lambda_min_ratio_transfer) ^ ((seq_len(nlam) - 1) / (nlam - 1))
    }
    
  }
  
  loss_trans=rep(0,times=length(lambda_beta))
  
  
  betaA_total=matrix(0,nrow=length(lambda_beta),ncol=p)
  
  for(i1 in 1:length(lambda_beta)){
    
    loss=rep(0,times=nfold)
    lambda_beta1=lambda_beta[i1]
    for(i2 in 1:nfold){
      
      c1=1:nrow(X_all)
      m1=round(nrow(X_all)/nfold,digits = 0)
      test=c1[(1+(i2-1)*m1):min((i2*m1),nrow(X_all))]
      X_all_train=X_all[-test,]
      y_all_train=y_all[-test]
      X_all_test=X_all[test,]
      y_all_test=y_all[test] 
      
      betaA=transfer_linear_fix(lambda = lambda_beta1,
                                X=X_all_train, y=y_all_train,
                                t0 = NULL,
                                beta0 = beta_start,
                                tol = tol_transfer,
                                maxit = maxit,
                                c = c_use,
                                r = 10,
                                use_linesearch = TRUE)
      betaA=betaA$beta
      loss[i2]=loss_linear(beta=betaA,y=y_all_test,X=X_all_test) 
    }
    loss_trans[i1]=mean(loss)
    
    betaA_total[i1,]=betaA
    
  }
  lambda_beta_use=max(lambda_beta[round(loss_trans,digits = 4)==min(round(loss_trans,digits = 4))])
  lambda_beta_index=max(c(1:length(lambda_beta))[round(loss_trans,digits = 4)==min(round(loss_trans,digits = 4))])
  betaA=transfer_linear_fix(lambda = lambda_beta_use,
                            X=X_all, y=y_all,
                            t0 = NULL,
                            beta0 = beta_start,
                            tol = tol_transfer,
                            maxit = maxit,
                            c = c_use,
                            r = 10,
                            use_linesearch = TRUE)
  beta_hat=betaA$beta
  #betaA=betaA_total[lambda_beta_index,]
  
  
  
  
  if(!intercept){
    return(list(beta_hat=c(0,beta_hat)))
  }
  if(intercept){
    return(list(beta_hat=beta_hat))
  }
  
  
}



#' Transfer learning for compositional linear regression under a zero-sum constraint
#'
#' The model supports non-compositional covariates plus
#' compositional covariates whose effects satisfy a zero-sum constraint.
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
CatlGLM_linear=function(target_data=NULL,source_data=NULL,nlam=100,Ncov=0,
                                     nfold=3,source_id="auto",intercept=TRUE,C=NULL,
                                     lambda_list=list(lambda_transfer=NULL,lambda_debias=NULL),maxit=600,
                                     tol_transfer=1e-6,tol_debias=1e-6,beta_start=NULL,delta_start=NULL,C0=NULL){
  .check_dataset(target_data, "target_data")
  .check_source_list(source_data, "source_data")
  .check_source_id(source_id)
  
  if(source_id=="auto"){

    LOSS1=matrix(0,nrow=nfold,ncol = length(source_data))
    
    LOSS1_lasso=numeric()
    beta_fold_transfer=list()
    for(i in 1:nfold){
      beta_fold_transfer[[i]]=list()
    }
    score_cv=matrix(0,nrow=nfold,ncol=length(source_data))
    beta_constrained_lasso=list()
    est_target=transfer_only_linear(target=target_data,
                                    lambda_beta=NULL,nfold=nfold,beta_start=NULL,maxit=400,tol_transfer=1e-5,Ncov=Ncov,nlam=50,intercept=intercept,C=C)  
    INFO=info_linear_n(y=target_data$y,X=target_data$x,beta=est_target$beta_hat,intercept=intercept)
    for(i in 1:nfold){
      
      fold_size=ceiling(length(target_data$y)/nfold)
      test=(1+(i-1)*fold_size):min(i*fold_size,length(target_data$y))
      train_x=target_data$x[-test,]
      train_y=target_data$y[-test]
      test_x=target_data$x[test,]
      test_y=target_data$y[test]
      #compute information matrix
      
      
      for(j in 1:length(source_data)){
        
        est_j=transfer_only_linear(target=list(x=rbind(train_x,source_data[[j]]$x),y=c(train_y,source_data[[j]]$y)),
                                     lambda_beta=NULL,nfold=nfold,beta_start=NULL,maxit=400,tol_transfer=1e-5,Ncov=Ncov,nlam=50,intercept=intercept,C=C)  
        beta_fold_transfer[[i]][[j]]=est_j$beta_hat
        LOSS1[i,j]=loss_linear(beta=est_j$beta_hat,X=cbind(1,test_x),y=test_y)
        
      }
      lambda_lasso=transfer_only_linear(target=list(x=train_x,y=train_y),lambda_beta=NULL,nfold=nfold,beta_start=NULL,maxit=400,tol_transfer=1e-5,Ncov=Ncov,nlam=50,intercept=intercept,C=C)  
      
      
      
      
      beta_lasso=lambda_lasso$beta_hat
      
      LOSS1_lasso[i]=loss_linear(y=test_y,X=cbind(1,test_x),beta=beta_lasso) 
      beta_constrained_lasso[[i]]=beta_lasso
      print(i)
    } 
    
    for(i in 1:nrow(score_cv)){
      for(j in 1:ncol(score_cv)){
        score_cv[i,j]=t(beta_fold_transfer[[i]][[j]]-beta_constrained_lasso[[i]])%*%INFO%*%(beta_fold_transfer[[i]][[j]]-beta_constrained_lasso[[i]])/LOSS1_lasso[i]
      }
    }
    LOSS_lasso=mean(LOSS1_lasso)
    LOSS_lasso_sd=sd(LOSS1_lasso)
    LOSS=colMeans(LOSS1)
    score_all=colMeans(score_cv)
    if(is.null(C0)){
      C0=min(score_all)
    }
    good_id=c(1:length(source_data))[LOSS-(1+C0)*LOSS_lasso<=0] 
    diff_source=LOSS-(1+C0)*LOSS_lasso
    good_source=source_data[good_id]
    if(length(good_source)==0){
      cat("no source are good\n", append = TRUE)   # write to .out file
      stop("Error: no valid IDs found length of 'good_id' is 0.")
      flush.console()               # ensure message appears immediately
    }
    final=transfer_linear(source=good_source,
                            target=target_data,
                            lambda_beta=lambda_list$lambda_transfer,lambda_delta=lambda_list$lambda_debias,nfold=nfold,beta_start=beta_start,
                            delta_start=delta_start,maxit=maxit,tol_transfer=tol_transfer,
                            tol_debias=tol_debias,Ncov=Ncov,nlam=nlam,C=C,intercept = intercept)
    
    
    return(list(beta_hat=final$beta_hat,transferrable_id=good_id,diff_source=diff_source))
  }
  if(source_id=="all"){
    
    final=transfer_linear(source=source_data,
                            target=target_data,
                            lambda_beta=lambda_list$lambda_transfer,lambda_delta=lambda_list$lambda_debias,nfold=nfold,beta_start=beta_start,
                            delta_start=delta_start,maxit=maxit,tol_transfer=tol_transfer,
                            tol_debias=tol_debias,Ncov=Ncov,nlam=nlam,C=C,intercept = intercept)
    return(list(beta_hat=final$beta_hat,transferrable_id=c(1:length(source_data))))
  }
  
  
}











