\name{CatlGLM_linear}
\alias{CatlGLM_linear}
\title{Transfer learning for compositional linear regression}
\description{
Fits the adaptive transfer-learning estimator for linear regression with
non-compositional covariates and compositional covariates subject to a zero-sum
constraint.
}
\usage{
CatlGLM_linear(
  target_data = NULL,
  source_data = NULL,
  nlam = 100,
  Ncov = 0,
  nfold = 3,
  source_id = "auto",
  intercept = TRUE,
  C = NULL,
  lambda_list = list(lambda_transfer = NULL, lambda_debias = NULL),
  maxit = 600,
  tol_transfer = 1e-6,
  tol_debias = 1e-6,
  beta_start = NULL,
  delta_start = NULL,
  C0 = NULL
)
}
\arguments{
  \item{target_data}{List with components \code{x} and \code{y}.}
  \item{source_data}{List of source datasets; each element must contain \code{x} and \code{y}.}
  \item{nlam}{Number of tuning parameters used in the transfer and debias steps.}
  \item{Ncov}{Number of non-compositional covariates.}
  \item{nfold}{Number of folds used in cross-validation.}
  \item{source_id}{Either \code{"auto"} or \code{"all"}.}
  \item{intercept}{Logical; whether to fit an intercept.}
  \item{C}{Optional linear-constraint vector.}
  \item{lambda_list}{List containing \code{lambda_transfer} and \code{lambda_debias}.}
  \item{maxit}{Maximum number of iterations.}
  \item{tol_transfer}{Convergence tolerance for the transfer step.}
  \item{tol_debias}{Convergence tolerance for the debias step.}
  \item{beta_start}{Optional starting vector for the transfer step.}
  \item{delta_start}{Optional starting vector for the debias step.}
  \item{C0}{Optional threshold for selecting transferable sources. If \code{NULL}, it is computed from the data.}
}
\value{
If \code{source_id = "all"}, returns a list with \code{beta_hat} and
\code{transferrable_id}. If \code{source_id = "auto"}, also returns
\code{diff_source}.
}
\examples{
set.seed(1)
nt <- 30
ns <- 30
p <- 4
xt <- matrix(rnorm(nt * p), nt, p)
xs <- matrix(rnorm(ns * p), ns, p)
yt <- xt[, 1] - xt[, 2] + rnorm(nt)
ys <- xs[, 1] - xs[, 2] + rnorm(ns)

target_data <- list(x = xt, y = yt)
source_data <- list(list(x = xs, y = ys))

# fit <- CatlGLM_linear(target_data, source_data, Ncov = 0, source_id = "all")
}
