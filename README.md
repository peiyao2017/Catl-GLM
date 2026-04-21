---
title: "CATL-GLM: a novel constraint-aware transfer learning method for compositional data under generalized linear model framework"
author: "Peiyao Wang, Weixi Chen, Jiyuan Hu"
---

# Overview

This R Markdown document reproduces the simulation setup and implementation of constraint-aware transfer learning for compositional data under generalized linear model framework (CATL-GLM). Users can run the code chunk by chunk or knit the full document after placing all required source files and `.RData` files in the working directory.

# Required files

Make sure the following files are available in your working directory:

 
- `effects_both_rare_and_abundant_binomial.RData`
- `effects_both_rare_and_abundant_linear.RData`

# Setup

Please install Rtools at https://cran.r-project.org/bin/windows/Rtools/ before installing the package.



```{r setup, message=FALSE, warning=FALSE}

library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(pROC)
library(coda4microbiome)
install.packages("remotes")
library(remotes)
remotes::install_github("JoJoWeixiChen1/Catl-GLM")
library(CatlGLM)
 

 
```

# 1. Binomial case

## Simulation parameters

```{r binomial-params}
rho <- 0.5 # inter-microbiome correlation
h1 <- 10 # transfer level
int <- 1 # intercept
have_intercept <- TRUE

n_test <- 100 # target test sample size
p <- 200 # number of compositional predictors
nlam <- 60 # number of tunning parameter
n_target <- 100 + n_test # total target sample size
n_source <- 100 # source sample size
constraint <- rep(1, times = p) # linear constraint vector
n_abund <- 0.1 * p # number of abundant microbiomes
```

## Mean and covariance structure of compositional features

```{r binomial-covariance}
mu_data <- c(rep(p / 2, times = n_abund), rep(2, times = p - n_abund))

sigmax <- matrix(0, nrow = p, ncol = p)
for (i in 1:nrow(sigmax)) {
  for (j in 1:ncol(sigmax)) {
    sigmax[i, j] <- rho^abs(i - j)
  }
}
```

## Load target effects

```{r binomial-load-effects}
load("C:/Users/wangp12/Downloads/make readme file/effects_both_rare_and_abundant_binomial.RData")
beta_target <- beta_total$beta_target
beta_target[1] <- int
```

## Generate source effects

The first `K0` sources are valid. The remaining sources are invalid.

```{r binomial-source-effects}
K <- 10 # total number of source data
K0 <- 5 # number of valid source

beta_source <- list()

for (i in 1:K0) {
  R1 <- rep(0, times = p / 2)
  heter_microb <- sample(1:length(R1), size = length(R1) / 5)
  R1_heter <- sample(c(1, -1), size = length(heter_microb), replace = TRUE, prob = c(0.5, 0.5))
  R1[heter_microb] <- R1_heter
  beta_source[[i]] <- beta_target + h1 / p * c(1, R1, -R1)
}

for (i in (K0 + 1):K) {
  microb_effect <- beta_target[-1]
  b <- rep(0, times = length(microb_effect) / 2)
  index <- (1:length(b))[microb_effect[1:(p / 2)] == 0 & microb_effect[(p / 2 + 1):p] == 0]
  s1 <- index
  s2 <- (1:length(b))[!(1:length(b)) %in% index]
  b[s1] <- runif(n = length(s1), min = 1.5, max = 2) +
    2 * h1 * sample(c(-1, 1), replace = TRUE, prob = c(0.5, 0.5), size = length(s1)) / (p / 2)
  b[s2] <- 2 + runif(n = length(s2), min = 1.5, max = 2) +
    2 * h1 * sample(c(-1, 1), replace = TRUE, prob = c(0.5, 0.5), size = length(s2)) / (p / 2)
  beta_source[[i]] <- c(1, b, -b)
}
```

## Generate target data

```{r binomial-target-data}
targetW <- matrix(0, nrow = ceiling(n_target), ncol = p) # microbial count matrix
targetX <- matrix(0, nrow = ceiling(n_target), ncol = p) # log relative abundance matrix
targetY <- rep(1, times = nrow(targetX)) # outcome vector

for (i in 1:nrow(targetX)) {
  b <- rep(0, times = ncol(targetX))
  a <- exp(mvrnorm(n = 1, mu = log(mu_data), Sigma = sigmax))

  for (j in 1:ncol(targetX)) {
    b[j] <- log(a[j] / sum(a))
  }

  prob_target <- exp(t(c(1, b)) %*% beta_target) / (1 + exp(t(c(1, b)) %*% beta_target))
  targetY[i] <- rbinom(n = 1, size = 1, prob = prob_target)
  targetW[i, ] <- a
  targetX[i, ] <- b
}
```

## Split target data into training and test sets

```{r binomial-split}
idx_test <- (nrow(targetX) - n_test + 1):nrow(targetX)
idx_train <- 1:(nrow(targetX) - n_test)

targetX_test <- targetX[idx_test, , drop = FALSE]
targetY_test <- targetY[idx_test]

targetX <- targetX[idx_train, , drop = FALSE]
targetY <- targetY[idx_train]
targetW <- targetW[idx_train, , drop = FALSE]
```

## Generate source data

```{r binomial-source-data}
sourceW <- vector("list", K)
sourceX <- vector("list", K)
sourceX_centered <- vector("list", K)
sourceY <- vector("list", K)

for (i in 1:K) {
  ep <- mvrnorm(n = 1, mu = rep(0, times = p), Sigma = 0.3 * diag(1, nrow = p, ncol = p))
  beta_source1 <- beta_source[[i]]
  SW <- matrix(0, nrow = floor(n_source), ncol = p)
  SX <- matrix(0, nrow = floor(n_source), ncol = p)
  SY <- rep(1, times = nrow(SX))

  for (i1 in 1:nrow(SX)) {
    b <- rep(0, times = ncol(SX))
    a <- exp(mvrnorm(n = 1, mu = log(mu_data), Sigma = sigmax + ep %*% t(ep)))

    for (j in 1:ncol(SX)) {
      b[j] <- log(a[j] / sum(a))
    }

    prob_source <- exp(t(c(1, b)) %*% beta_source1) / (1 + exp(t(c(1, b)) %*% beta_source1))
    SY[i1] <- rbinom(n = 1, size = 1, prob = prob_source)
    SW[i1, ] <- a
    SX[i1, ] <- b
  }

  sourceW[[i]] <- SW
  sourceX[[i]] <- SX
  sourceY[[i]] <- SY
}
```

## Organize target and source data

```{r binomial-assemble}
target_total <- list(y = targetY, x = targetX)

source_total <- vector("list", length(sourceX))
for (i in 1:length(sourceX)) {
  source_total[[i]] <- list(x = sourceX[[i]], y = sourceY[[i]])
}
```

## Run CATL when valid sources are known

```{r binomial-est-true}
est_true <- CatlGLM_binomial(
  target_data = target_total,
  source_data = source_total[1:K0],
  nlam = 100, Ncov = 0, nfold = 3, source_id = "all",
  intercept = have_intercept, C0 = 0.5, C = NULL,
  lambda_list = list(lambda_transfer = NULL, lambda_debias = NULL),
  maxit = 500, tol_transfer = 1e-5, tol_debias = 1e-5,
  beta_start = NULL, delta_start = NULL
)
```

## Run CATL when valid sources are unknown

This step may take a long time.

```{r binomial-est-detect}
est_detect <- CatlGLM_binomial(
  target_data = target_total,
  source_data = source_total,
  nlam = 100, Ncov = 0, nfold = 3, source_id = "auto",
  intercept = have_intercept, C0 = NULL, C = NULL,
  lambda_list = list(lambda_transfer = NULL, lambda_debias = NULL),
  maxit = 500, tol_transfer = 1e-5, tol_debias = 1e-5,
  beta_start = NULL, delta_start = NULL
)
```

## Extract point estimates

```{r binomial-coef}
beta_est_true <- est_true$beta_hat
beta_est_detect <- est_detect$beta_hat
```

## Construct confidence intervals using selected transferable sources

```{r binomial-ci}
inf_est <- CatlGLM_inf(
  family = "binomial", intercept = FALSE, target = target_total,
  source = source_total[est_detect$transferrable_id],
  nodewise.transfer.source.id = "all", level = 0.95,
  beta.hat = est_detect$beta_hat
)
```

## Compute AUC on the test set

```{r binomial-auc}
prob_est_true <- exp(cbind(1, targetX_test) %*% beta_est_true) / (1 + exp(cbind(1, targetX_test) %*% beta_est_true))
prob_est_detect <- exp(cbind(1, targetX_test) %*% beta_est_detect) / (1 + exp(cbind(1, targetX_test) %*% beta_est_detect))

auc_est_true <- roc(response = targetY_test, predictor = prob_est_true[, 1])$auc
auc_est_detect <- roc(response = targetY_test, predictor = prob_est_detect[, 1])$auc

auc_est_true
auc_est_detect
```

# 2. Linear case

## Simulation parameters

```{r linear-params}
rho <- 0.5 # inter-microbiome correlation
h1 <- 10 # transfer level
int <- 1 # intercept
have_intercept <- FALSE
sd_err <- 6 # dispersion parameter

n_test <- 100 # target test sample size
p <- 200 # number of compositional predictors
nlam <- 60 # number of tunning parameter
n_target <- 100 + n_test # total target sample size
n_source <- 100 # source sample size
constraint <- rep(1, times = p) # linear constraint vector
n_abund <- 0.1 * p # number of abundant microbiomes
```

## Mean and covariance structure

```{r linear-covariance}
mu_data <- c(rep(p / 2, times = n_abund), rep(2, times = p - n_abund))

sigmax <- matrix(0, nrow = p, ncol = p)
for (i in 1:nrow(sigmax)) {
  for (j in 1:ncol(sigmax)) {
    sigmax[i, j] <- rho^abs(i - j)
  }
}
```

## Load target effects

```{r linear-load-effects}
load("C:/Users/wangp12/Downloads/make readme file/effects_both_rare_and_abundant_linear.RData")
beta_target <- beta_total$beta_target
beta_target[1] <- int
```

## Generate source effects

```{r linear-source-effects}
K <- 40 # total number of source data
K0 <- 5 # number of valid source

beta_source <- list()

for (i in 1:K0) {
  R1 <- rep(0, times = p / 2)
  heter_microb <- sample(1:length(R1), size = length(R1) / 5)
  R1_heter <- sample(c(1, -1), size = length(heter_microb), replace = TRUE, prob = c(0.5, 0.5))
  R1[heter_microb] <- R1_heter
  beta_source[[i]] <- beta_target + h1 / p * c(1, R1, -R1)
}

for (i in (K0 + 1):K) {
  microb_effect <- beta_target[-1]
  b <- rep(0, times = length(microb_effect) / 2)
  index <- (1:length(b))[microb_effect[1:(p / 2)] == 0 & microb_effect[(p / 2 + 1):p] == 0]
  s1 <- index
  s2 <- (1:length(b))[!(1:length(b)) %in% index]
  b[s1] <- runif(n = length(s1), min = 1.5, max = 2) +
    2 * h1 * sample(c(-1, 1), replace = TRUE, prob = c(0.5, 0.5), size = length(s1)) / (p / 2)
  b[s2] <- 2 + runif(n = length(s2), min = 1.5, max = 2) +
    2 * h1 * sample(c(-1, 1), replace = TRUE, prob = c(0.5, 0.5), size = length(s2)) / (p / 2)
  beta_source[[i]] <- c(1, b, -b)
}
```

## Generate target data

```{r linear-target-data}
targetW <- matrix(0, nrow = ceiling(n_target), ncol = p) # microbial count matrix
targetX <- matrix(0, nrow = ceiling(n_target), ncol = p) # log relative abundance matrix

for (i in 1:nrow(targetX)) {
  b <- rep(0, times = ncol(targetX))
  a <- exp(mvrnorm(n = 1, mu = log(mu_data), Sigma = sigmax))

  for (j in 1:ncol(targetX)) {
    b[j] <- log(a[j] / sum(a))
  }

  targetW[i, ] <- a
  targetX[i, ] <- b
}

targetY <- cbind(rep(1, nrow(targetX)), targetX) %*% beta_target +
  rnorm(n = n_target, mean = 0, sd = sd_err)
```

## Split target data into training and test sets

```{r linear-split}
idx_test <- (nrow(targetX) - n_test + 1):nrow(targetX)
idx_train <- 1:(nrow(targetX) - n_test)

targetX_test <- targetX[idx_test, , drop = FALSE]
targetY_test <- targetY[idx_test]

targetX <- targetX[idx_train, , drop = FALSE]
targetY <- targetY[idx_train]
targetW <- targetW[idx_train, , drop = FALSE]
```

## Generate source data

```{r linear-source-data}
sourceW <- vector("list", K)
sourceX <- vector("list", K)
sourceX_centered <- vector("list", K)
sourceY <- vector("list", K)

for (i in 1:K) {
  ep <- mvrnorm(n = 1, mu = rep(0, times = p), Sigma = 0.3 * diag(1, nrow = p, ncol = p))
  beta_source1 <- beta_source[[i]]
  SW <- matrix(0, nrow = floor(n_source), ncol = p)
  SX <- matrix(0, nrow = floor(n_source), ncol = p)

  for (i1 in 1:nrow(SX)) {
    b <- rep(0, times = ncol(SX))
    a <- exp(mvrnorm(n = 1, mu = log(mu_data), Sigma = sigmax + ep %*% t(ep)))

    for (j in 1:ncol(SX)) {
      b[j] <- log(a[j] / sum(a))
    }

    SW[i1, ] <- a
    SX[i1, ] <- b
  }

  SY <- cbind(1, SX) %*% beta_source[[i]] + rnorm(n = n_source, mean = 0, sd = sd_err)

  sourceW[[i]] <- SW
  sourceX[[i]] <- SX
  sourceY[[i]] <- SY
}
```

## Organize target and source data

```{r linear-assemble}
target_total <- list(y = targetY, x = targetX)

source_total <- vector("list", length(sourceX))
for (i in 1:length(sourceX)) {
  source_total[[i]] <- list(x = sourceX[[i]], y = sourceY[[i]])
}
```

## Run transfer learning when valid sources are known

```{r linear-est-true}
est_true <- CatlGLM_linear(
  target_data = target_total,
  source_data = source_total[1:K0],
  nlam = 100, Ncov = 0, nfold = 3, source_id = "all",
  intercept = have_intercept, C0 = 0.5, C = NULL,
  lambda_list = list(lambda_transfer = NULL, lambda_debias = NULL),
  maxit = 500, tol_transfer = 1e-5, tol_debias = 1e-5,
  beta_start = NULL, delta_start = NULL
)
```

## Run transfer learning when valid sources are unknown

This step may take a long time.

```{r linear-est-detect}
est_detect <- CatlGLM_linear(
  target_data = target_total,
  source_data = source_total,
  nlam = 100, Ncov = 0, nfold = 3, source_id = "auto",
  intercept = have_intercept, C0 = NULL, C = NULL,
  lambda_list = list(lambda_transfer = NULL, lambda_debias = NULL),
  maxit = 500, tol_transfer = 1e-5, tol_debias = 1e-5,
  beta_start = NULL, delta_start = NULL
)
```

## Extract point estimates

```{r linear-coef}
beta_est_true <- est_true$beta_hat
beta_est_detect <- est_detect$beta_hat
```

## Estimate residual scale and standardize data for inference

```{r linear-scale}
sigma_est <- sqrt(mean((targetY - cbind(1, targetX) %*% beta_est_detect)^2))

target_total_scaled_est <- list(
  y = target_total$y / sigma_est,
  x = target_total$x / sigma_est
)

source_total_scaled_est <- list()
for (i in 1:length(source_total)) {
  source_total_scaled_est[[i]] <- list(
    y = source_total[[i]]$y / sigma_est,
    x = source_total[[i]]$x / sigma_est
  )
}
```

## Construct confidence intervals using selected transferable sources

```{r linear-ci}
inf_est <- CatlGLM_inf(
  family = "gaussian", intercept = FALSE, target = target_total_scaled_est,
  source = source_total_scaled_est[est_detect$transferrable_id],
  nodewise.transfer.source.id = "all", level = 0.95,
  beta.hat = est_detect$beta_hat
)
```

## Compute prediction error on the test set

```{r linear-mse}
mse_est_detect <- sqrt(
  loss_linear(
    beta = beta_est_detect,
    y = targetY_test,
    X = cbind(rep(1, nrow(targetX_test)), targetX_test)
  )
)

mse_est_detect
```

# Session info

```{r session-info}
sessionInfo()
```
