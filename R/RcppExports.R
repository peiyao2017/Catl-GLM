# Generated wrapper functions for compiled code.

g_ls_cpp <- function(beta, X, y) {
  .Call(`_CatlGLM_g_ls_cpp`, beta, X, y)
}

grad_g_ls_cpp <- function(beta, X, y) {
  .Call(`_CatlGLM_grad_g_ls_cpp`, beta, X, y)
}

prox_l1_zerosum_rank1_cpp <- function(v, tau, c) {
  .Call(`_CatlGLM_prox_l1_zerosum_rank1_cpp`, v, tau, c)
}

G_t_ls_constrained_cpp <- function(yvec, t, X, y, lambda, c) {
  .Call(`_CatlGLM_G_t_ls_constrained_cpp`, yvec, t, X, y, lambda, c)
}

prox_l1_shift_zerosum_rank1_cpp <- function(v, tau, c, betaA) {
  .Call(`_CatlGLM_prox_l1_shift_zerosum_rank1_cpp`, v, tau, c, betaA)
}

G_t_ls_debias_constrained_cpp <- function(yvec, t, X, y, lambda, c, betaA) {
  .Call(`_CatlGLM_G_t_ls_debias_constrained_cpp`, yvec, t, X, y, lambda, c, betaA)
}

proj_ker_ct_cpp <- function(u, c) {
  .Call(`_CatlGLM_proj_ker_ct_cpp`, u, c)
}

soft_thresh_cpp <- function(v, tau) {
  .Call(`_CatlGLM_soft_thresh_cpp`, v, tau)
}

g_cpp <- function(beta, X, y) {
  .Call(`_CatlGLM_g_cpp`, beta, X, y)
}

grad_g_cpp <- function(beta, X, y) {
  .Call(`_CatlGLM_grad_g_cpp`, beta, X, y)
}

soft_thresh_debias_cpp <- function(v, betaA, tau) {
  .Call(`_CatlGLM_soft_thresh_debias_cpp`, v, betaA, tau)
}

G_t_cpp <- function(yvec, t, X, y, lambda, c) {
  .Call(`_CatlGLM_G_t_cpp`, yvec, t, X, y, lambda, c)
}

G_t_debias_cpp <- function(yvec, t, X, y, lambda, c, betaA) {
  .Call(`_CatlGLM_G_t_debias_cpp`, yvec, t, X, y, lambda, c, betaA)
}
