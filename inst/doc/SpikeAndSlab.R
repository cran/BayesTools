## ----setup, include = FALSE---------------------------------------------------
# is_check <- ("CheckExEnv" %in% search()) ||
#              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv()))
is_check <- F
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  eval     = !is_check,
  dev      = "png"
)
if(.Platform$OS.type == "windows"){
  knitr::opts_chunk$set(dev.args = list(type = "cairo"))
}

## -----------------------------------------------------------------------------
library(BayesTools)

## -----------------------------------------------------------------------------
set.seed(-68) # set seed for reproducibility

N     <- 100      # number of observations
x     <- rnorm(N) # continuous predictor
alpha <- -0.5     # intercept
beta  <- 0.15     # small effect

# compute the mean parameter for each predictor value
mu <- alpha + beta * x

# generate the response for each observation 
y  <- rnorm(N, mean = mu, sd = 1) 

## -----------------------------------------------------------------------------
summary(lm(y ~ x))

## -----------------------------------------------------------------------------
model_likelihood <- 
"model{
  for(i in 1:N){
    y[i] ~ dnorm(mu[i], pow(sigma, -2))
  }
}
"

## -----------------------------------------------------------------------------
formula_M0 <- list("mu" = ~ 1)
formula_M1 <- list("mu" = ~ 1 + x)

## -----------------------------------------------------------------------------
# prior on the test parameter
prior_beta  <- prior(distribution = "normal", parameters = list(mean = 0, sd = 0.5))

# priors on the nuisance parameters
prior_int   <- prior(distribution = "normal", parameters = list(mean = 0, sd = 5))
prior_sigma <- prior(distribution = "normal", parameters = list(mean = 0, sd = 5), truncation = list(0, Inf))

# the data list
data_list <- list(
  y = y,
  N = N
)
data_formula <- data.frame(
  x = x
)

## -----------------------------------------------------------------------------
M0 <- JAGS_fit(
  # specification for the `model_likelihood` part
  model_syntax = model_likelihood,
  data         = list(y = y, N = N),
  prior_list   = list("sigma" = prior_sigma),

  # specification for the `formula_M0` part 
  formula_list       = formula_M0,
  formula_prior_list = list("mu" = list("intercept" = prior_int)),
  formula_data_list  = list("mu" = data_formula),
  
  # seed for reproducibility
  seed         = 0
)

M1 <- JAGS_fit(
  model_syntax = model_likelihood,
  data         = list(y = y, N = N),
  prior_list   = list("sigma" = prior_sigma),
  formula_list       = formula_M1,
  formula_prior_list = list("mu" = list("intercept" = prior_int, "x" = prior_beta)),
  formula_data_list  = list("mu" = data_formula),
  seed         = 1
)

## -----------------------------------------------------------------------------
JAGS_estimates_table(M1)

## -----------------------------------------------------------------------------
log_posterior <- function(parameters, data){
  sum(dnorm(
    x    = data[["y"]],
    mean = parameters[["mu"]],
    sd   = parameters[["sigma"]],
    log  = TRUE
  ))
}

## -----------------------------------------------------------------------------
marglik_model_H0 <- JAGS_bridgesampling(
  # specification for the model part
  fit           = M0,
  log_posterior = log_posterior,
  data          = list(y = y, N = N),
  prior_list    = list("sigma" = prior_sigma),

  # specification for the formula` part 
  formula_list       = formula_M0,
  formula_prior_list = list("mu" = list("intercept" = prior_int)),
  formula_data_list  = list("mu" = data_formula)
)

marglik_model_H1 <- JAGS_bridgesampling(
  fit           = M1,
  log_posterior = log_posterior,
  data          = list(y = y, N = N),
  prior_list    = list("sigma" = prior_sigma),
  formula_list       = formula_M1,
  formula_prior_list = list("mu" = list("intercept" = prior_int, "x" = prior_beta)),
  formula_data_list  = list("mu" = data_formula),
)

## -----------------------------------------------------------------------------
models_list <- models_inference(list(
  list(model = M0, marglik = marglik_model_H0, prior_weights = 1/2),
  list(model = M1, marglik = marglik_model_H1, prior_weights = 1/2)
))
ensemble_info <- ensemble_inference(models_list, parameters = "x", is_null_list = list("x" = c(TRUE, FALSE)))

ensemble_inference_table(ensemble_info, parameters = "x")

## -----------------------------------------------------------------------------
prior_beta_spike_and_slab <- prior_spike_and_slab(
  prior_parameter = prior(distribution = "normal", parameters = list(mean = 0, sd = 0.5)),
  prior_inclusion = prior(distribution = "spike", parameters = list(location = 0.5)) 
)


## -----------------------------------------------------------------------------
MS <- JAGS_fit(
  model_syntax = model_likelihood,
  data         = list(y = y, N = N),
  prior_list   = list("sigma" = prior_sigma),
  formula_list       = formula_M1,
  formula_prior_list = list("mu" = list("intercept" = prior_int, "x" = prior_beta_spike_and_slab)),
  formula_data_list  = list("mu" = data_formula),
  seed         = 1
)

## -----------------------------------------------------------------------------
JAGS_estimates_table(MS, conditional = TRUE)

## -----------------------------------------------------------------------------
JAGS_inference_table(MS)

