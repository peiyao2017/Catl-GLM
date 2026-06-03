---
title: "CATL-GLM: a novel constraint-aware transfer learning method for compositional data under generalized linear model framework"
author: "Peiyao Wang, Weixi Chen, Jiyuan Hu"
---

# Overview

This R Markdown document reproduces the simulation setup and implementation of constraint-aware transfer learning for compositional data under generalized linear model framework (CATL-GLM). Users can run the code chunk by chunk or knit the full document after placing all required source files and `.RData` files in the working directory. The introduction of data harmonization is in section 3.

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
# Install from GitHub
remotes::install_github("peiyao2017/Catl-GLM")
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

# 3. Data harmonization

This section gives a reproducible example of data harmonization for transfer learning with microbiome data and continuous outcome. To perform data harmonization using the CatlGLM package, all datasets must follow the same rule of variable names, and all missing values should be represented by "NA". The workflow includes:

1. Generating artificial target and source microbiome datasets.
2. Aligning target and source datasets by common variable names.
3. Removing observations with missing outcomes.
4. Processing microbiome features into log-relative abundance.
5. Removing variables and samples with excessive missingness.
6. Imputing missing non-microbiome covariates using `missForest`.
7. Returning the processed data in a transfer-learning-ready format.

The final output contains a target list and a source list. Each dataset is represented by two elements:

- `y`: the outcome vector.
- `x`: the predictor data frame containing covariates and log-relative-abundance microbiome features.

## Required package

```{r load-package}
library(missForest)
```

## Generate artificial microbiome data

Artificial microbiome and non-microbiome covariates are generated for one target dataset and ten source datasets. Each dataset has mostly common variables, but also some non-common variables. Missing values are represented by `NA`, and zero microbiome relative abundances are also introduced.

```{r generate-data}
set.seed(123)

K <- 10                # number of source datasets
p <- 50                # number of abundant microbiomes
q <- 50                # number of rare microbiomes
M <- p + q

n_target <- 100        # target sample size
n_source <- rep(100, K) # source sample sizes

zero_rate <- 0.05       # proportion of zero microbiome counts
cov_missing_rate <- 0.05 # missing rate of covariates
ra_missing_rate <- 0.05  # missing rate of relative abundance values

# Common and candidate microbiome names
common_microbes <- paste0("microb", 1:75)
microbe_pool <- paste0("microb", 1:150)

# Covariate names
common_covariates <- c("cov1", "cov2")
all_covariates <- paste0("cov", 1:5)

# Store feature names for target and source datasets
microbe_names_list <- vector("list", K + 1)
cov_names_list <- vector("list", K + 1)

for (k in 1:(K + 1)) {
  n_microbes_k <- sample(seq(80, 110, by = 2), 1)

  extra_microbes <- sample(
    setdiff(microbe_pool, common_microbes),
    size = n_microbes_k - length(common_microbes)
  )

  microbe_names_list[[k]] <- c(common_microbes, extra_microbes)

  n_cov_k <- sample(3:5, 1)

  extra_covariates <- sample(
    setdiff(all_covariates, common_covariates),
    size = n_cov_k - length(common_covariates)
  )

  cov_names_list[[k]] <- c(common_covariates, extra_covariates)
}

# Generate global target beta
beta_cov_global <- rnorm(5, mean = 0, sd = 3)
names(beta_cov_global) <- all_covariates

b <- rnorm(M / 2, mean = 0, sd = 1)
beta_microbe_global <- c(-b, b)
names(beta_microbe_global) <- paste0("microb", 1:M)

beta_global <- c(beta_cov_global, beta_microbe_global)

# Containers
all_count_data <- vector("list", K + 1)
all_logRA_data <- vector("list", K + 1)
all_final_data <- vector("list", K + 1)
all_beta <- vector("list", K + 1)

for (k in 1:(K + 1)) {

  if (k == 1) {
    n <- n_target
  } else {
    n <- n_source[k - 1]
  }

  microbe_names <- microbe_names_list[[k]]
  cov_names <- cov_names_list[[k]]

  n_microbes <- length(microbe_names)
  n_cov <- length(cov_names)
  m <- n_microbes + n_cov

  # Generate microbiome log-counts
  microbe_id <- as.integer(sub("microb", "", microbe_names))

  mu_vec <- ifelse(
    microbe_id <= p,
    20,
    1
  )

  log_count_mat <- matrix(
    rnorm(n * n_microbes, mean = rep(mu_vec, each = n), sd = 1),
    nrow = n,
    ncol = n_microbes
  )

  count_mat <- exp(log_count_mat)
  colnames(count_mat) <- microbe_names

  # Generate covariates
  cov_dat <- data.frame(matrix(NA, nrow = n, ncol = n_cov))
  colnames(cov_dat) <- cov_names

  for (cov in cov_names) {
    if (cov %in% c("cov1", "cov2")) {
      cov_dat[[cov]] <- rbinom(n, size = 1, prob = 0.5)
    } else {
      cov_dat[[cov]] <- rnorm(n, mean = 0, sd = 1)
    }
  }

  # Auxiliary data: log-relative abundance from original counts
  RA_clean <- count_mat / rowSums(count_mat)
  logRA_clean <- log(RA_clean)

  # Construct X for outcome generation
  X_clean <- cbind(as.matrix(cov_dat), logRA_clean)

  # Target beta restricted to current features
  beta_cov <- beta_cov_global[cov_names]

  beta_microbe <- beta_microbe_global[microbe_names]
  beta_microbe[is.na(beta_microbe)] <- 0

  beta_current <- c(beta_cov, beta_microbe)
  names(beta_current) <- colnames(X_clean)

  # Source beta perturbation
  if (k > 1) {

    if (n_microbes %% 2 != 0) {
      stop("Number of microbiome features must be even for zero-sum source perturbation.")
    }

    r_vec <- sample(
      c(-1, 1),
      size = n_microbes / 2,
      replace = TRUE,
      prob = c(0.5, 0.5)
    )

    R_microbe <- c(-r_vec, r_vec)

    beta_current[microbe_names] <- beta_current[microbe_names] +
      40 / m * R_microbe
  }

  all_beta[[k]] <- beta_current

  # Generate outcome
  epsilon <- rnorm(n, mean = 0, sd = 6)
  Y <- as.vector(X_clean %*% beta_current + epsilon)

  # Randomly make some counts zero
  zero_index <- matrix(
    rbinom(n * n_microbes, size = 1, prob = zero_rate),
    nrow = n,
    ncol = n_microbes
  )

  count_mat_zero <- count_mat
  count_mat_zero[zero_index == 1] <- 0

  # Convert zero-added counts to relative abundance
  RA_obs <- count_mat_zero / rowSums(count_mat_zero)
  RA_obs[is.nan(RA_obs)] <- 0

  # Randomly make some relative abundance values NA
  RA_missing_index <- matrix(
    rbinom(n * n_microbes, size = 1, prob = ra_missing_rate),
    nrow = n,
    ncol = n_microbes
  )

  RA_obs[RA_missing_index == 1] <- NA

  # Randomly make some covariate values NA
  cov_obs <- cov_dat

  for (j in 1:ncol(cov_obs)) {
    miss_index <- rbinom(n, size = 1, prob = cov_missing_rate)
    cov_obs[miss_index == 1, j] <- NA
  }

  # Final observed dataset: outcome, covariates, and relative abundance
  final_dat <- data.frame(
    Y = Y,
    cov_obs,
    RA_obs,
    check.names = FALSE
  )

  all_count_data[[k]] <- count_mat_zero
  all_logRA_data[[k]] <- logRA_clean
  all_final_data[[k]] <- final_dat
}

target_data <- all_final_data[[1]]

source_data <- all_final_data[-1]
names(source_data) <- paste0("source", 1:K)

beta_target <- all_beta[[1]]
beta_source <- all_beta[-1]
names(beta_source) <- paste0("source", 1:K)
```

## Initial alignment of target and source datasets

The function `align_TL_data()` conducts initial data processing. It removes observations with missing outcome values, keeps only common features across all datasets, and aligns the target and source datasets by variable name.

```{r align-function}
align_TL_data <- function(target, source, outcome) {

  if (!is.data.frame(target)) {
    stop("target must be a data.frame.")
  }

  if (!is.list(source) || !all(sapply(source, is.data.frame))) {
    stop("source must be a list of data.frames.")
  }

  all_data <- c(list(target = target), source)

  # Check outcome exists
  for (i in seq_along(all_data)) {
    if (!outcome %in% colnames(all_data[[i]])) {
      stop(paste0(
        "Outcome column '", outcome,
        "' is missing in dataset ", i, "."
      ))
    }
  }

  # Remove observations with missing outcomes
  removed_outcome_missing <- vector("list", length(all_data))

  for (i in seq_along(all_data)) {
    dat <- all_data[[i]]

    removed_outcome_missing[[i]] <- which(is.na(dat[[outcome]]))
    dat <- dat[!is.na(dat[[outcome]]), , drop = FALSE]

    all_data[[i]] <- dat
  }

  # Move outcome to first column
  all_data <- lapply(all_data, function(dat) {
    dat[, c(outcome, setdiff(colnames(dat), outcome)), drop = FALSE]
  })

  # Find common features excluding outcome
  feature_sets <- lapply(all_data, function(dat) {
    setdiff(colnames(dat), outcome)
  })

  common_features <- Reduce(intersect, feature_sets)

  # Use target feature order
  target_feature_order <- setdiff(
    colnames(all_data[[1]]),
    outcome
  )

  common_features <- target_feature_order[
    target_feature_order %in% common_features
  ]

  # Align all datasets
  all_data <- lapply(all_data, function(dat) {
    dat[, c(outcome, common_features), drop = FALSE]
  })

  target_aligned <- all_data[[1]]
  source_aligned <- all_data[-1]

  if (is.null(names(source)) || any(names(source) == "")) {
    names(source_aligned) <- paste0("source", seq_along(source_aligned))
  } else {
    names(source_aligned) <- names(source)
  }

  return(list(
    target = target_aligned,
    source = source_aligned,
    common_features = common_features,
    removed_outcome_missing = removed_outcome_missing
  ))
}
```
#Initial alignment
The function align_TL_data conduct initial data processing: only common variables are retained for each data, each data is aligned by variable name and observations with missing outcome values are removed.

```{r run-alignment}
aligned_data <- align_TL_data(
  target = target_data,
  source = source_data,
  outcome = "Y"
)

target_aligned <- aligned_data$target #aligned target data
source_aligned <- aligned_data$source #aligned source data
```

## Microbiome data processing and covariate imputation

The function `process_microbiome_TL_data()` performs the final harmonization step. It first removes samples with too many missing predictor values, then removes covariates and microbiome features with too much missingness. Since each dataset may lose different variables, the function then retains only common covariates and common microbiome features across all datasets.

After removing variables with too many missings, Microbiome variables are processed according to the argument `microbiome.scale`:

- If `microbiome.scale = "count"` or `"relative_abundance"`, remaining microbiome `NA` values are replaced by zero, then all zeros are replaced by a small positive pseudo-count, values are renormalized to unit row sum, and log-relative abundance is calculated. Missing non-microbiome covariates are imputed by `missForest`, using both covariates and log-relative-abundance microbiome features as predictors.
 
- If `microbiome.scale = "log_relative_abundance"`, values are first converted back to relative abundance by exponentiation, then above processing is performed.



```{r process-function}
process_microbiome_TL_data <- function(
  target,
  source,
  outcome,
  covariate,
  microbiome.scale = c("count", "relative_abundance", "log_relative_abundance"),
  microbe_missing_rate = 0.2,
  covariate_missing_rate = 0.5,
  sample_missing_rate = 0.5,
  pseudo.count = NULL,
  use.missForest = TRUE
) {

  microbiome.scale <- match.arg(microbiome.scale)

  if (use.missForest && !requireNamespace("missForest", quietly = TRUE)) {
    stop("Please install missForest first: install.packages('missForest')")
  }

  all_data <- c(list(target = target), source)

  for (k in seq_along(all_data)) {
    if (!outcome %in% colnames(all_data[[k]])) {
      stop(paste0("Outcome column is missing in dataset ", k, "."))
    }
  }

  processed_each <- vector("list", length(all_data))
  dropped_samples <- vector("list", length(all_data))
  dropped_covariates <- vector("list", length(all_data))
  dropped_microbes <- vector("list", length(all_data))

  for (k in seq_along(all_data)) {
    dat <- all_data[[k]]

    cov_k <- intersect(covariate, colnames(dat))
    microbe_k <- setdiff(colnames(dat), c(outcome, covariate))

    dat <- dat[, c(outcome, cov_k, microbe_k), drop = FALSE]

    # 1. Remove samples with too many missing predictor values
    miss_mat <- is.na(dat[, c(cov_k, microbe_k), drop = FALSE])
    sample_miss <- rowMeans(miss_mat)
    keep_sample <- sample_miss <= sample_missing_rate

    dropped_samples[[k]] <- which(!keep_sample)
    dat <- dat[keep_sample, , drop = FALSE]

    # 2. Remove covariates with too much missingness
    if (length(cov_k) > 0) {
      cov_miss <- colMeans(is.na(dat[, cov_k, drop = FALSE]))
      keep_cov <- cov_miss <= covariate_missing_rate
      dropped_covariates[[k]] <- names(cov_miss)[!keep_cov]
      cov_k <- names(cov_miss)[keep_cov]
    } else {
      dropped_covariates[[k]] <- character(0)
    }

    # 3. Remove microbiomes with too much missingness
    if (length(microbe_k) > 0) {
      microbe_miss <- colMeans(is.na(dat[, microbe_k, drop = FALSE]))
      keep_microbe <- microbe_miss <= microbe_missing_rate
      dropped_microbes[[k]] <- names(microbe_miss)[!keep_microbe]
      microbe_k <- names(microbe_miss)[keep_microbe]
    } else {
      dropped_microbes[[k]] <- character(0)
    }

    processed_each[[k]] <- dat[, c(outcome, cov_k, microbe_k), drop = FALSE]
  }

  # 4. Retain common covariates and common microbiomes
  cov_sets <- lapply(processed_each, function(dat) {
    intersect(covariate, colnames(dat))
  })

  microbe_sets <- lapply(processed_each, function(dat) {
    setdiff(colnames(dat), c(outcome, covariate))
  })

  common_cov <- Reduce(intersect, cov_sets)
  common_microbe <- Reduce(intersect, microbe_sets)

  common_cov <- covariate[covariate %in% common_cov]

  target_microbe_order <- setdiff(
    colnames(processed_each[[1]]),
    c(outcome, covariate)
  )

  common_microbe <- target_microbe_order[
    target_microbe_order %in% common_microbe
  ]

  processed_each <- lapply(processed_each, function(dat) {
    dat[, c(outcome, common_cov, common_microbe), drop = FALSE]
  })

  dataset_id <- rep(seq_along(processed_each), sapply(processed_each, nrow))
  combined <- do.call(rbind, processed_each)

  Y <- combined[[outcome]]
  cov_dat <- combined[, common_cov, drop = FALSE]
  microb_dat <- as.matrix(combined[, common_microbe, drop = FALSE])
  storage.mode(microb_dat) <- "numeric"

  # 5. Convert microbiome data to relative abundance scale
  if (microbiome.scale == "count") {
    microb_dat[is.na(microb_dat)] <- 0
    microb_dat[microb_dat < 0] <- 0

    row_sum <- rowSums(microb_dat)
    keep_nonzero <- row_sum > 0

    Y <- Y[keep_nonzero]
    cov_dat <- cov_dat[keep_nonzero, , drop = FALSE]
    microb_dat <- microb_dat[keep_nonzero, , drop = FALSE]
    dataset_id <- dataset_id[keep_nonzero]

    RA <- microb_dat / rowSums(microb_dat)
  }

  if (microbiome.scale == "relative_abundance") {
    microb_dat[is.na(microb_dat)] <- 0
    microb_dat[microb_dat < 0] <- 0

    row_sum <- rowSums(microb_dat)
    keep_nonzero <- row_sum > 0

    Y <- Y[keep_nonzero]
    cov_dat <- cov_dat[keep_nonzero, , drop = FALSE]
    microb_dat <- microb_dat[keep_nonzero, , drop = FALSE]
    dataset_id <- dataset_id[keep_nonzero]

    RA <- microb_dat / rowSums(microb_dat)
  }

  if (microbiome.scale == "log_relative_abundance") {
    RA_raw <- exp(microb_dat)
    RA_raw[is.na(RA_raw)] <- 0
    RA_raw[RA_raw < 0] <- 0

    row_sum <- rowSums(RA_raw)
    keep_nonzero <- row_sum > 0

    Y <- Y[keep_nonzero]
    cov_dat <- cov_dat[keep_nonzero, , drop = FALSE]
    RA_raw <- RA_raw[keep_nonzero, , drop = FALSE]
    dataset_id <- dataset_id[keep_nonzero]

    RA <- RA_raw / rowSums(RA_raw)
  }

  # 6. Replace zero RA by pseudo-count, renormalize, then log-transform
  if (is.null(pseudo.count)) {
    pseudo.count <- max(
      min(RA[RA > 0], na.rm = TRUE) / 2,
      1e-10
    )
  }

  RA[RA == 0] <- pseudo.count
  RA <- RA / rowSums(RA)

  logRA_df <- as.data.frame(log(RA))
  colnames(logRA_df) <- common_microbe

  # 7. Impute covariates only, using covariates + logRA as predictors
  cov_imp <- cov_dat

  if (use.missForest && ncol(cov_imp) > 0 && any(is.na(cov_imp))) {
    X_for_impute <- data.frame(
      cov_imp,
      logRA_df,
      check.names = FALSE
    )

    for (v in common_cov) {
      obs_vals <- na.omit(unique(X_for_impute[[v]]))

      if (length(obs_vals) <= 2 && all(obs_vals %in% c(0, 1))) {
        X_for_impute[[v]] <- factor(X_for_impute[[v]], levels = c(0, 1))
      }
    }

    X_imp <- missForest::missForest(X_for_impute)$ximp

    cov_imp <- X_imp[, common_cov, drop = FALSE]
  }

  # 8. Reconstruct processed datasets
  combined_processed <- data.frame(
    cov_imp,
    logRA_df,
    check.names = FALSE
  )

  processed_split <- split(
    data.frame(
      dataset_id = dataset_id,
      Y = Y,
      combined_processed,
      check.names = FALSE
    ),
    dataset_id
  )

  target_df <- processed_split[[as.character(1)]]

  target_data <- list(
    y = target_df[[outcome]],
    x = target_df[, setdiff(colnames(target_df), c("dataset_id", outcome)), drop = FALSE]
  )

  source_data <- vector("list", length(source))

  for (k in seq_along(source)) {
    source_df <- processed_split[[as.character(k + 1)]]

    source_data[[k]] <- list(
      y = source_df[[outcome]],
      x = source_df[, setdiff(colnames(source_df), c("dataset_id", outcome)), drop = FALSE]
    )
  }

  if (is.null(names(source)) || any(names(source) == "")) {
    names(source_data) <- paste0("source", seq_along(source_data))
  } else {
    names(source_data) <- names(source)
  }

  return(list(
    target_data = target_data,
    source_data = source_data,
    outcome = outcome,
    covariates = common_cov,
    microbiome_features = common_microbe,
    dropped_samples = dropped_samples,
    dropped_covariates = dropped_covariates,
    dropped_microbes = dropped_microbes,
    pseudo.count = pseudo.count
  ))
}
```

## Apply the processing function

The following code applies the microbiome processing function to the aligned target and source datasets. The output format matches the structure commonly required by transfer-learning functions.

```{r apply-processing}
harmonized_data <- process_microbiome_TL_data(
  target = target_aligned,
  source = source_aligned,
  outcome = "Y",
  covariate = c("cov1", "cov2"),
  microbiome.scale = "relative_abundance",
  microbe_missing_rate = 0.2,
  covariate_missing_rate = 0.5,
  sample_missing_rate = 0.5,
  pseudo.count = NULL,
  use.missForest = TRUE
)

target_data_harmonized <- harmonized_data$target_data
source_data_harmonized <- harmonized_data$source_data
```

## Inspect the output

The target dataset has an outcome vector `y` and a predictor data frame `x`.

```{r inspect-target}
str(target_data_harmonized)
dim(target_data_harmonized$x)
length(target_data_harmonized$y)
```

Each source dataset has the same structure.

```{r inspect-source}
names(source_data_harmonized)
str(source_data_harmonized$source1)
dim(source_data_harmonized$source1$x)
length(source_data_harmonized$source1$y)
```

The following objects record which variables or samples were removed during processing.

```{r inspect-removed}
harmonized_data$dropped_samples
harmonized_data$dropped_covariates
harmonized_data$dropped_microbes
harmonized_data$pseudo.count
```

## Notes

For supervised transfer learning, observations with missing outcomes are removed rather than imputed. Missing microbiome values are handled on the relative-abundance scale before log transformation, while missing non-microbiome covariates are imputed using both covariates and microbiome log-relative abundance as predictors.



## Session info

```{r session-info}
sessionInfo()
```
