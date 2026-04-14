set.seed(1)

nt <- 50
ns1 <- 40
ns2 <- 40
p <- 5

xt <- matrix(rnorm(nt * p), nt, p)
xs1 <- matrix(rnorm(ns1 * p), ns1, p)
xs2 <- matrix(rnorm(ns2 * p), ns2, p)

# Linear example
beta_true <- c(1, -1, 0.5, -0.5, 0)
yt_lin <- xt %*% beta_true + rnorm(nt)
ys1_lin <- xs1 %*% beta_true + rnorm(ns1)
ys2_lin <- xs2 %*% beta_true + rnorm(ns2)

lin_target <- list(x = xt, y = as.numeric(yt_lin))
lin_source <- list(
  list(x = xs1, y = as.numeric(ys1_lin)),
  list(x = xs2, y = as.numeric(ys2_lin))
)

fit_lin <- CatlGLM_linear(lin_target, lin_source, source_id = "all")

ci_lin <- CatlGLM_inf(
  target = lin_target,
  source = lin_source,
  family = "gaussian",
  beta.hat = fit_lin$beta_hat,
  cores = 1,
  level = 0.95,
  intercept = TRUE,
  ncov = 0
)

ci_lin$CI
ci_lin$b.hat
ci_lin$var.est

## Linear example, source = auto
set.seed(2)

nt <- 80
ns1 <- 60
ns2 <- 60
p <- 5

xt  <- matrix(rnorm(nt * p),  nt,  p)
xs1 <- matrix(rnorm(ns1 * p), ns1, p)
xs2 <- matrix(rnorm(ns2 * p), ns2, p)

beta_true <- c(1, -1, 0.5, -0.5, 0)

beta_bad  <- c(-1, 1, -0.5, 0.5, 1)

# ------------------------
# Linear example
# ------------------------
yt_lin  <- xt  %*% beta_true + rnorm(nt,  sd = 1)
ys1_lin <- xs1 %*% beta_true + rnorm(ns1, sd = 1)
ys2_lin <- xs2 %*% beta_bad  + rnorm(ns2, sd = 1.5)

lin_target <- list(x = xt, y = as.numeric(yt_lin))
lin_source <- list(
  list(x = xs1, y = as.numeric(ys1_lin)),  
  list(x = xs2, y = as.numeric(ys2_lin))   
)

fit_lin_auto <- CatlGLM_linear(
  lin_target,
  lin_source,
  source_id = "auto",
  intercept = TRUE
)




# Logistic example
eta_t <- xt %*% beta_true
eta_s1 <- xs1 %*% beta_true
eta_s2 <- xs2 %*% beta_true

yt_bin <- rbinom(nt, 1, plogis(eta_t))
ys1_bin <- rbinom(ns1, 1, plogis(eta_s1))
ys2_bin <- rbinom(ns2, 1, plogis(eta_s2))

bin_target <- list(x = xt, y = yt_bin)
bin_source <- list(
  list(x = xs1, y = ys1_bin),
  list(x = xs2, y = ys2_bin)
)

fit_bin <- CatlGLM_binomial(bin_target, bin_source, source_id = "all")

ci_bin <- CatlGLM_inf(
  target = bin_target,
  source = bin_source,
  family = "binomial",
  beta.hat = fit_bin$beta_hat,
  cores = 1,
  level = 0.95,
  intercept = TRUE,
  ncov = 0
)


ci_bin$CI
ci_bin$b.hat
ci_bin$var.est



## Logistic example, source = auto
eta_t  <- xt  %*% beta_true
eta_s1 <- xs1 %*% beta_true
eta_s2 <- xs2 %*% beta_bad

yt_bin  <- rbinom(nt,  1, plogis(eta_t))
ys1_bin <- rbinom(ns1, 1, plogis(eta_s1))
ys2_bin <- rbinom(ns2, 1, plogis(eta_s2))

bin_target <- list(x = xt, y = yt_bin)
bin_source <- list(
  list(x = xs1, y = ys1_bin),  
  list(x = xs2, y = ys2_bin)  
)

fit_bin_auto <- CatlGLM_binomial(
  bin_target,
  bin_source,
  source_id = "auto",
  intercept = TRUE
)


