\name{CatlGLM_inf}
\alias{CatlGLM_inf}
\title{Confidence intervals for compositional transfer-learning estimates}
\description{
Computes de-biased point estimates and Wald confidence intervals for the
estimators returned by \code{CatlGLM_linear()} and
\code{CatlGLM_binomial()}.
}
\usage{
CatlGLM_inf(
  target = NULL,
  source = NULL,
  family = c("gaussian", "binomial"),
  beta.hat = NULL,
  nodewise.transfer.source.id = "all",
  cores = 1,
  level = 0.95,
  intercept = TRUE,
  ncov = 0
)
}
\arguments{
  \item{target}{Target dataset as a list with components \code{x} and \code{y}.}
  \item{source}{List of source datasets; each element must contain \code{x} and \code{y}.}
  \item{family}{Either \code{"gaussian"} or \code{"binomial"}.}
  \item{beta.hat}{Coefficient vector returned by \code{CatlGLM_linear()} or \code{CatlGLM_binomial()}.}
  \item{nodewise.transfer.source.id}{Source indices used in nodewise transfer; default is \code{"all"}.}
  \item{cores}{Number of cores used for parallel computation.}
  \item{level}{Confidence level.}
  \item{intercept}{Logical; must match the fitted model.}
  \item{ncov}{Number of non-compositional covariates.}
}
\value{
A list containing \code{b.hat}, \code{beta.hat}, \code{CI}, and
\code{var.est}.
}
\examples{
# Typically used after fitting CatlGLM_linear() or CatlGLM_binomial().
}
