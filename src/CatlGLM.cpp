
#include <Rcpp.h>
using namespace Rcpp;

// =========================
// Linear-regression helpers
// =========================
double g_ls_cpp(NumericVector beta, NumericMatrix X, NumericVector y) {
  int n = X.nrow(), p = X.ncol();
  double out = 0.0;
  for (int i = 0; i < n; i++) {
    double xb = 0.0;
    for (int j = 0; j < p; j++) xb += X(i, j) * beta[j];
    double r = xb - y[i];
    out += r * r;
  }
  return 0.5 * out / n;
}

NumericVector grad_g_ls_cpp(NumericVector beta, NumericMatrix X, NumericVector y) {
  int n = X.nrow(), p = X.ncol();
  NumericVector r(n);
  for (int i = 0; i < n; i++) {
    double xb = 0.0;
    for (int j = 0; j < p; j++) xb += X(i, j) * beta[j];
    r[i] = xb - y[i];
  }
  NumericVector grad(p);
  for (int j = 0; j < p; j++) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += X(i, j) * r[i];
    grad[j] = s / n;
  }
  return grad;
}

inline double soft_scalar(double x, double tau) {
  if (x > tau) return x - tau;
  if (x < -tau) return x + tau;
  return 0.0;
}

NumericVector prox_l1_zerosum_rank1_cpp(NumericVector v, double tau, NumericVector c) {
  int p = v.size();
  if (c.size() != p) stop("v and c must have same length");

  std::vector<int> idx;
  idx.reserve(p);
  for (int j = 0; j < p; ++j) if (c[j] != 0.0) idx.push_back(j);
  int m = static_cast<int>(idx.size());

  if (m == 0) {
    NumericVector out(p);
    for (int j = 0; j < p; ++j) out[j] = soft_scalar(v[j], tau);
    return out;
  }

  auto f = [&](double mu) -> double {
    double s = 0.0;
    for (int k = 0; k < m; ++k) {
      int j = idx[k];
      s += soft_scalar(v[j] - mu, tau);
    }
    return s;
  };

  double lo = -1.0, hi = 1.0;
  double flo = f(lo), fhi = f(hi);
  int it = 0;
  while (flo * fhi > 0.0 && it < 60) {
    lo *= 2.0;
    hi *= 2.0;
    ++it;
    flo = f(lo);
    fhi = f(hi);
  }

  double mu = 0.0;
  if (flo * fhi <= 0.0) {
    for (int t = 0; t < 80; ++t) {
      mu = 0.5 * (lo + hi);
      double fm = f(mu);
      if (fm == 0.0) break;
      if (flo * fm < 0.0) {
        hi = mu;
        fhi = fm;
      } else {
        lo = mu;
        flo = fm;
      }
    }
  }

  NumericVector beta(p);
  for (int j = 0; j < p; ++j) {
    if (c[j] != 0.0) beta[j] = soft_scalar(v[j] - mu, tau);
    else beta[j] = soft_scalar(v[j], tau);
  }
  return beta;
}

NumericVector G_t_ls_constrained_cpp(NumericVector yvec, double t,
                                     NumericMatrix X, NumericVector y,
                                     double lambda, NumericVector c) {
  int p = yvec.size();
  NumericVector gy = grad_g_ls_cpp(yvec, X, y);
  NumericVector v(p);
  for (int j = 0; j < p; ++j) v[j] = yvec[j] - t * gy[j];
  NumericVector beta_t = prox_l1_zerosum_rank1_cpp(v, t * lambda, c);
  NumericVector G(p);
  for (int j = 0; j < p; ++j) G[j] = (yvec[j] - beta_t[j]) / t;
  return G;
}

inline double soft0(double x, double t) {
  if (x > t) return x - t;
  if (x < -t) return x + t;
  return 0.0;
}

NumericVector prox_l1_shift_zerosum_rank1_cpp(NumericVector v,
                                              double tau,
                                              NumericVector c,
                                              NumericVector betaA) {
  int p = v.size();
  if (c.size() != p || betaA.size() != p) stop("v, c, betaA must have same length");

  NumericVector w(p);
  double target = 0.0;
  std::vector<int> idx;
  idx.reserve(p);
  for (int j = 0; j < p; ++j) {
    w[j] = v[j] - betaA[j];
    if (c[j] != 0.0) {
      idx.push_back(j);
      target -= betaA[j];
    }
  }
  int m = static_cast<int>(idx.size());

  if (m == 0) {
    NumericVector z(p), beta(p);
    for (int j = 0; j < p; ++j) z[j] = soft0(w[j], tau);
    for (int j = 0; j < p; ++j) beta[j] = betaA[j] + z[j];
    return beta;
  }

  auto f = [&](double mu) -> double {
    double s = 0.0;
    for (int k = 0; k < m; ++k) {
      int j = idx[k];
      s += soft0(w[j] - mu, tau);
    }
    return s - target;
  };

  double lo = -1.0, hi = 1.0;
  double flo = f(lo), fhi = f(hi);
  int it = 0;
  while (flo * fhi > 0.0 && it < 60) {
    lo *= 2.0;
    hi *= 2.0;
    ++it;
    flo = f(lo);
    fhi = f(hi);
  }

  double mu = 0.0;
  if (flo * fhi <= 0.0) {
    for (int t = 0; t < 80; ++t) {
      mu = 0.5 * (lo + hi);
      double fm = f(mu);
      if (fm == 0.0) break;
      if (flo * fm < 0.0) {
        hi = mu;
        fhi = fm;
      } else {
        lo = mu;
        flo = fm;
      }
    }
  }

  NumericVector z(p), beta(p);
  for (int j = 0; j < p; ++j) {
    if (c[j] != 0.0) z[j] = soft0(w[j] - mu, tau);
    else z[j] = soft0(w[j], tau);
    beta[j] = betaA[j] + z[j];
  }
  return beta;
}

NumericVector G_t_ls_debias_constrained_cpp(NumericVector yvec, double t,
                                            NumericMatrix X, NumericVector y,
                                            double lambda, NumericVector c,
                                            NumericVector betaA) {
  int p = yvec.size();
  NumericVector gy = grad_g_ls_cpp(yvec, X, y);
  NumericVector v(p);
  for (int j = 0; j < p; ++j) v[j] = yvec[j] - t * gy[j];
  NumericVector beta_t = prox_l1_shift_zerosum_rank1_cpp(v, t * lambda, c, betaA);
  NumericVector G(p);
  for (int j = 0; j < p; ++j) G[j] = (yvec[j] - beta_t[j]) / t;
  return G;
}

// ===========================
// Logistic-regression helpers
// ===========================
NumericVector proj_ker_ct_cpp(NumericVector u, NumericVector c) {
  int p = u.size();
  if (c.size() != p) stop("u and c must have the same length");
  double denom = 0.0, num = 0.0;
  for (int i = 0; i < p; i++) {
    denom += c[i] * c[i];
    num += c[i] * u[i];
  }
  if (denom <= 1e-15) stop("c must not be all zeros");
  double alpha = num / denom;
  NumericVector out(p);
  for (int i = 0; i < p; i++) out[i] = u[i] - alpha * c[i];
  return out;
}

NumericVector soft_thresh_cpp(NumericVector v, double tau) {
  int p = v.size();
  NumericVector out(p);
  for (int i = 0; i < p; i++) {
    double x = v[i];
    double ax = std::abs(x);
    out[i] = (ax > tau) ? ((x > 0 ? 1 : -1) * (ax - tau)) : 0.0;
  }
  return out;
}

double g_cpp(NumericVector beta, NumericMatrix X, NumericVector y) {
  int n = X.nrow(), p = X.ncol();
  double out = 0.0;
  for (int i = 0; i < n; i++) {
    double eta = 0.0;
    for (int j = 0; j < p; j++) eta += X(i, j) * beta[j];
    if (eta > 40) eta = 40;
    if (eta < -40) eta = -40;
    out += log1p(std::exp(eta)) - y[i] * eta;
  }
  return out / n;
}

NumericVector grad_g_cpp(NumericVector beta, NumericMatrix X, NumericVector y) {
  int n = X.nrow(), p = X.ncol();
  NumericVector grad(p);
  for (int i = 0; i < n; i++) {
    double eta = 0.0;
    for (int j = 0; j < p; j++) eta += X(i, j) * beta[j];
    if (eta > 40) eta = 40;
    if (eta < -40) eta = -40;
    double p_i = 1.0 / (1.0 + std::exp(-eta));
    double diff = p_i - y[i];
    for (int j = 0; j < p; j++) grad[j] += X(i, j) * diff;
  }
  for (int j = 0; j < p; j++) grad[j] /= n;
  return grad;
}

NumericVector soft_thresh_debias_cpp(NumericVector v, NumericVector betaA, double tau) {
  int p = v.size();
  if (betaA.size() != p) stop("v and betaA must have the same length");
  NumericVector out(p);
  for (int j = 0; j < p; j++) {
    double vj = v[j];
    double bAj = betaA[j];
    if (vj - tau > bAj) {
      out[j] = vj - tau;
    } else if (vj + tau < bAj) {
      out[j] = vj + tau;
    } else {
      out[j] = bAj;
    }
  }
  return out;
}

NumericVector G_t_cpp(NumericVector yvec, double t,
                      NumericMatrix X, NumericVector y,
                      double lambda, NumericVector c) {
  NumericVector gy = grad_g_cpp(yvec, X, y);
  int p = yvec.size();
  NumericVector v(p);
  for (int i = 0; i < p; i++) v[i] = yvec[i] - t * gy[i];
  NumericVector u = soft_thresh_cpp(v, t * lambda);
  NumericVector beta_t = proj_ker_ct_cpp(u, c);
  NumericVector out(p);
  for (int i = 0; i < p; i++) out[i] = (yvec[i] - beta_t[i]) / t;
  return out;
}

NumericVector G_t_debias_cpp(NumericVector yvec,
                             double t,
                             NumericMatrix X,
                             NumericVector y,
                             double lambda,
                             NumericVector c,
                             NumericVector betaA) {
  NumericVector gy = grad_g_cpp(yvec, X, y);
  int p = yvec.size();
  NumericVector v(p);
  for (int i = 0; i < p; i++) v[i] = yvec[i] - t * gy[i];
  NumericVector u = soft_thresh_debias_cpp(v, betaA, t * lambda);
  NumericVector beta_t = proj_ker_ct_cpp(u, c);
  NumericVector out(p);
  for (int i = 0; i < p; i++) out[i] = (yvec[i] - beta_t[i]) / t;
  return out;
}

// =========================
// .Call wrappers + register
// =========================
#define CatlGLM_BEGIN_RCPP Rcpp::RObject rcpp_result_gen; Rcpp::RNGScope rcpp_rngScope_gen;
#define CatlGLM_END_RCPP return rcpp_result_gen;

extern "C" SEXP _CatlGLM_g_ls_cpp(SEXP betaSEXP, SEXP XSEXP, SEXP ySEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector beta(betaSEXP);
  NumericMatrix X(XSEXP);
  NumericVector y(ySEXP);
  rcpp_result_gen = wrap(g_ls_cpp(beta, X, y));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_grad_g_ls_cpp(SEXP betaSEXP, SEXP XSEXP, SEXP ySEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector beta(betaSEXP);
  NumericMatrix X(XSEXP);
  NumericVector y(ySEXP);
  rcpp_result_gen = wrap(grad_g_ls_cpp(beta, X, y));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_prox_l1_zerosum_rank1_cpp(SEXP vSEXP, SEXP tauSEXP, SEXP cSEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector v(vSEXP);
  double tau = as<double>(tauSEXP);
  NumericVector c(cSEXP);
  rcpp_result_gen = wrap(prox_l1_zerosum_rank1_cpp(v, tau, c));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_G_t_ls_constrained_cpp(SEXP yvecSEXP, SEXP tSEXP, SEXP XSEXP, SEXP ySEXP, SEXP lambdaSEXP, SEXP cSEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector yvec(yvecSEXP);
  double t = as<double>(tSEXP);
  NumericMatrix X(XSEXP);
  NumericVector y(ySEXP);
  double lambda = as<double>(lambdaSEXP);
  NumericVector c(cSEXP);
  rcpp_result_gen = wrap(G_t_ls_constrained_cpp(yvec, t, X, y, lambda, c));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_prox_l1_shift_zerosum_rank1_cpp(SEXP vSEXP, SEXP tauSEXP, SEXP cSEXP, SEXP betaASEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector v(vSEXP);
  double tau = as<double>(tauSEXP);
  NumericVector c(cSEXP);
  NumericVector betaA(betaASEXP);
  rcpp_result_gen = wrap(prox_l1_shift_zerosum_rank1_cpp(v, tau, c, betaA));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_G_t_ls_debias_constrained_cpp(SEXP yvecSEXP, SEXP tSEXP, SEXP XSEXP, SEXP ySEXP, SEXP lambdaSEXP, SEXP cSEXP, SEXP betaASEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector yvec(yvecSEXP);
  double t = as<double>(tSEXP);
  NumericMatrix X(XSEXP);
  NumericVector y(ySEXP);
  double lambda = as<double>(lambdaSEXP);
  NumericVector c(cSEXP);
  NumericVector betaA(betaASEXP);
  rcpp_result_gen = wrap(G_t_ls_debias_constrained_cpp(yvec, t, X, y, lambda, c, betaA));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_proj_ker_ct_cpp(SEXP uSEXP, SEXP cSEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector u(uSEXP);
  NumericVector c(cSEXP);
  rcpp_result_gen = wrap(proj_ker_ct_cpp(u, c));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_soft_thresh_cpp(SEXP vSEXP, SEXP tauSEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector v(vSEXP);
  double tau = as<double>(tauSEXP);
  rcpp_result_gen = wrap(soft_thresh_cpp(v, tau));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_g_cpp(SEXP betaSEXP, SEXP XSEXP, SEXP ySEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector beta(betaSEXP);
  NumericMatrix X(XSEXP);
  NumericVector y(ySEXP);
  rcpp_result_gen = wrap(g_cpp(beta, X, y));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_grad_g_cpp(SEXP betaSEXP, SEXP XSEXP, SEXP ySEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector beta(betaSEXP);
  NumericMatrix X(XSEXP);
  NumericVector y(ySEXP);
  rcpp_result_gen = wrap(grad_g_cpp(beta, X, y));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_soft_thresh_debias_cpp(SEXP vSEXP, SEXP betaASEXP, SEXP tauSEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector v(vSEXP);
  NumericVector betaA(betaASEXP);
  double tau = as<double>(tauSEXP);
  rcpp_result_gen = wrap(soft_thresh_debias_cpp(v, betaA, tau));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_G_t_cpp(SEXP yvecSEXP, SEXP tSEXP, SEXP XSEXP, SEXP ySEXP, SEXP lambdaSEXP, SEXP cSEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector yvec(yvecSEXP);
  double t = as<double>(tSEXP);
  NumericMatrix X(XSEXP);
  NumericVector y(ySEXP);
  double lambda = as<double>(lambdaSEXP);
  NumericVector c(cSEXP);
  rcpp_result_gen = wrap(G_t_cpp(yvec, t, X, y, lambda, c));
  CatlGLM_END_RCPP
}

extern "C" SEXP _CatlGLM_G_t_debias_cpp(SEXP yvecSEXP, SEXP tSEXP, SEXP XSEXP, SEXP ySEXP, SEXP lambdaSEXP, SEXP cSEXP, SEXP betaASEXP) {
  CatlGLM_BEGIN_RCPP
  NumericVector yvec(yvecSEXP);
  double t = as<double>(tSEXP);
  NumericMatrix X(XSEXP);
  NumericVector y(ySEXP);
  double lambda = as<double>(lambdaSEXP);
  NumericVector c(cSEXP);
  NumericVector betaA(betaASEXP);
  rcpp_result_gen = wrap(G_t_debias_cpp(yvec, t, X, y, lambda, c, betaA));
  CatlGLM_END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
  {"_CatlGLM_g_ls_cpp", (DL_FUNC) &_CatlGLM_g_ls_cpp, 3},
  {"_CatlGLM_grad_g_ls_cpp", (DL_FUNC) &_CatlGLM_grad_g_ls_cpp, 3},
  {"_CatlGLM_prox_l1_zerosum_rank1_cpp", (DL_FUNC) &_CatlGLM_prox_l1_zerosum_rank1_cpp, 3},
  {"_CatlGLM_G_t_ls_constrained_cpp", (DL_FUNC) &_CatlGLM_G_t_ls_constrained_cpp, 6},
  {"_CatlGLM_prox_l1_shift_zerosum_rank1_cpp", (DL_FUNC) &_CatlGLM_prox_l1_shift_zerosum_rank1_cpp, 4},
  {"_CatlGLM_G_t_ls_debias_constrained_cpp", (DL_FUNC) &_CatlGLM_G_t_ls_debias_constrained_cpp, 7},
  {"_CatlGLM_proj_ker_ct_cpp", (DL_FUNC) &_CatlGLM_proj_ker_ct_cpp, 2},
  {"_CatlGLM_soft_thresh_cpp", (DL_FUNC) &_CatlGLM_soft_thresh_cpp, 2},
  {"_CatlGLM_g_cpp", (DL_FUNC) &_CatlGLM_g_cpp, 3},
  {"_CatlGLM_grad_g_cpp", (DL_FUNC) &_CatlGLM_grad_g_cpp, 3},
  {"_CatlGLM_soft_thresh_debias_cpp", (DL_FUNC) &_CatlGLM_soft_thresh_debias_cpp, 3},
  {"_CatlGLM_G_t_cpp", (DL_FUNC) &_CatlGLM_G_t_cpp, 6},
  {"_CatlGLM_G_t_debias_cpp", (DL_FUNC) &_CatlGLM_G_t_debias_cpp, 7},
  {NULL, NULL, 0}
};

extern "C" void R_init_CatlGLM(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
