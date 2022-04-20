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
data("kitchen_rolls")
x <- kitchen_rolls$mean_NEO[kitchen_rolls$rotation == "counter"]
y <- kitchen_rolls$mean_NEO[kitchen_rolls$rotation == "clock"]

## ----fig_dist, fig.width = 4, fig.height = 4, dpi = 300, out.width = "60%", fig.align = "center"----
h1 <- hist(x, breaks = 15, plot = FALSE)
h2 <- hist(y, breaks = 15, plot = FALSE)
par(mar = c(4, 4, 0, 1))
plot(h1, col = rgb(0,0,1,1/4), xlim = c(-1, 2), ylim = c(0, 16), las = 1, main = "", xlab = "mean NEO PI-R")
plot(h2, col = rgb(1,0,0,1/4), add = TRUE)
legend("topright", legend = c("Counter", "Clock"), fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), bty = "n")

## -----------------------------------------------------------------------------
ttest_model <- 
"model{
  for(i in 1:Nx){
    x[i] ~ dnorm(mu - delta*sigma/2, pow(sigma, -2))
  }
  for(i in 1:Ny){
    y[i] ~ dnorm(mu + delta*sigma/2, pow(sigma, -2))
  }
}
"

## -----------------------------------------------------------------------------
ttest_priors_H0 <- list(
  mu    = prior("cauchy",      parameters = list(location = 0, scale = 10)),
  sigma = prior("exponential", parameters = list(rate = 2)),
  delta = prior("spike",       parameters = list(location = 0))
)
ttest_priors_Hp <- list(
  mu    = prior("cauchy",      parameters = list(location = 0, scale = 10)),
  sigma = prior("exponential", parameters = list(rate = 2)),
  delta = prior("cauchy",      parameters = list(location = 0, scale = 1),
                truncation = list(lower = 0, upper = Inf))
)

## ----fig_priors, fig.width = 8, fig.height = 8, dpi = 300, out.width = "100%", fig.align = "center"----
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(ttest_priors_H0$mu,    par_name = bquote(mu), xlim = c(-50, 50))
plot(ttest_priors_H0$sigma, par_name = bquote(sigma), xlim = c(0, 50))
plot(ttest_priors_H0$delta, par_name = bquote(delta), 
     xlim = c(-1, 1), main = bquote(H[0]))
plot(ttest_priors_Hp$delta, par_name = bquote(delta), 
     xlim = c(0, 5), main = bquote(H[1]))

## -----------------------------------------------------------------------------
ttest_data <- list(
  x  = x,
  y  = y,
  Nx = length(x),
  Ny = length(y)
)

## -----------------------------------------------------------------------------
ttest_model_H0 <- JAGS_fit(
  model_syntax = ttest_model,
  data         = ttest_data,
  prior_list   = ttest_priors_H0,
  seed         = 0
) 

ttest_model_Hp <- JAGS_fit(
  model_syntax = ttest_model,
  data         = ttest_data,
  prior_list   = ttest_priors_Hp,
  seed         = 1
) 

## -----------------------------------------------------------------------------
runjags_estimates_table(ttest_model_Hp)

## -----------------------------------------------------------------------------
log_posterior <- function(parameters, data){
  loglik_x <- sum(dnorm(
    x    = data[["x"]],
    mean = parameters[["mu"]] - parameters[["delta"]] * parameters[["sigma"]] / 2,
    sd   = parameters[["sigma"]],
    log  = TRUE
  ))
  loglik_y <- sum(dnorm(
    x    = data[["y"]],
    mean = parameters[["mu"]] + parameters[["delta"]] * parameters[["sigma"]] / 2,
    sd   = parameters[["sigma"]], 
    log  = TRUE
  ))
  return(loglik_x + loglik_y)
}

## -----------------------------------------------------------------------------
marglik_model_H0 <- JAGS_bridgesampling(
  fit           = ttest_model_H0,
  log_posterior = log_posterior,
  data          = ttest_data,
  prior_list    = ttest_priors_H0
)

marglik_model_Hp <- JAGS_bridgesampling(
  fit           = ttest_model_Hp,
  log_posterior = log_posterior,
  data          = ttest_data,
  prior_list    = ttest_priors_Hp
)

## -----------------------------------------------------------------------------
bridgesampling::bf(marglik_model_H0, marglik_model_Hp)

## -----------------------------------------------------------------------------
models_list <- models_inference(list(
  list(model = ttest_model_H0, marglik = marglik_model_H0, prior_weights = 1/2),
  list(model = ttest_model_Hp, marglik = marglik_model_Hp, prior_weights = 1/2)
))
enseble_info <- ensemble_inference(models_list, parameters = "delta", is_null_list = list("delta" = c(TRUE, FALSE)))

ensemble_inference_table(enseble_info, "delta", BF01 = TRUE)

## -----------------------------------------------------------------------------
BayesFactor_ttest <- BayesFactor::ttestBF(x = x, y = y, rscale = "wide", nullInterval = c(0, Inf))
BayesFactor_ttest
1/exp(BayesFactor_ttest@bayesFactor$bf[2])

