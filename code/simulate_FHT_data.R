library(statmod)

run_nlm <- function(X, Z, times, observed) {
  # X and Z matrices
  d <- dim(X)[2]
  p <- dim(Z)[2]
  N <- dim(X)[1]
  X_design_matrix <- cbind(rep(1, N), X)
  Z_design_matrix <- cbind(rep(1, N), Z)
  censored_survival_times <- times

  loglikfunc <- function(params) {
    beta_size <- d + 1
    gamma_size <- p + 1
    beta_ <- params[1:beta_size]
    gamma_ <- params[(beta_size+1):(beta_size+gamma_size)]
    # standardize X
    y0 <- exp(X_design_matrix %*% beta_)
    mu <- Z_design_matrix %*% gamma_
    times <- censored_survival_times
    log_f <- log(y0) - 0.5*log(2*pi*times^3) - ((y0 + mu*times)^2)/(2*times)
    log_S <- log(pnorm((y0+mu*times)/sqrt(times)) - exp(-2*y0*mu)*pnorm((mu*times-y0)/sqrt(times)))
    loglik <- sum(observed*log_f) + sum((1 - observed)*log_S)
    total_loglikelihood <- sum(loglik)
    return(total_loglikelihood)
  }
  
  minus_loglik <- function(params) {
    value <- -loglikfunc(params)
    return(value)
  }
  
  N_params <- p+d+2
  initial_params <- runif(N_params, min=1, max=2)
  print(minus_loglik(initial_params))
  fit <- nlm(minus_loglik, initial_params)
  fit
}
### SIMULATION EXPERIMENT
set.seed(2)
N <- 1000
d <- 1
p <- 1
beta_ <- c(2, -1)
gamma_ <- c(-2, -1)
X <- cbind(c(rep(1, 500), rep(0, 500)))
Z <- cbind(c(rep(1, 200), rep(4, 300), rep(10, 500)))
y_intercepts <- rep(1, N)
X_design_matrix <- cbind(y_intercepts, X)
Z_design_matrix <- cbind(y_intercepts, Z)
y0 <- exp(X_design_matrix %*% beta_)
mu <- Z_design_matrix %*% gamma_
sigma_2 <- 1
mu_IG <- - y0/mu
lambda_IG <- (y0/sigma_2)^2
survival_times <- rinvgauss(N, mu_IG, lambda_IG)
W <- 0.7 # censoring time
censored_survival_times <- survival_times
censored_survival_times[is.na(survival_times)] <- W # censor!
observed <- ifelse(censored_survival_times < W, 1, 0)
times <- pmin(survival_times, W)

## Done simulating, now estimate

fit <- run_nlm(X, Z, censored_survival_times, observed)
beta_hat <- fit$estimate[1:2]
gamma_hat <- fit$estimate[3:4]
beta_
beta_hat
gamma_
gamma_hat







kaplanmeiers <- function(t_vector, observed) {
  sorted_t_vector <- times[order(t_vector)]
  sorted_observed <- observed[order(t_vector)]
  at_risk <- calculate_at_risk(sorted_t_vector)
  N <- length(t_vector)
  out <- rep(1, N)
  out[2:N] <- out[1:(N-1)] * (1 - sorted_observed[1:(N-1)]/at_risk[1:(N-1)])
  return(out)
}



S_FHT <- function(t, mu, y0) {
  pnorm((y0+mu*t)/sqrt(t)) - exp(-2*y0*mu)*pnorm((mu*t-y0)/sqrt(t))
}

getminusFHT_loglik <- function(t, delta) {
  FHT_loglik <- function(params) {
    y0 <- params[1]
    mu <- params[2]
    log_f <- log(y0) - 0.5*log(2*pi*t^3) - (y0 + mu*t)^2/(2*t)
    log_S <- log(pnorm((mu*t + y0)/sqrt(t)) - exp(-2*y0*mu)*pnorm((mu*t - y0))/sqrt(t))
    sum(delta*log_f + sum(1-delta)*log_S)
  }
  minusFHT_loglik <- function(params) {
    -FHT_loglik(params)
  }
  minusFHT_loglik
}

set.seed(1)
N <- 100
y0 <- 4
mu <- -1
mu_IG <- - y0/mu
lambda_IG <- y0^2
survival_times <- rinvgauss(N, mu_IG, lambda_IG)
W <- 6
censored_survival_times <- survival_times
censored_survival_times[is.na(survival_times)] <- W # censor!
observed <- ifelse(censored_survival_times < W, 1, 0)
times <- pmin(censored_survival_times, W)
sorted_times <- sort(times)
S_hats <- kaplanmeiers(times, observed)
plot(sorted_times, S_hats, xlab='days', typ='s')


minusFHT_loglik <- getFHT_loglik(sorted_times, sorted_observed)
params <- nlm(FHT_loglik, c(0.2, 0.2))
y0_hat <- params$estimate[1]
mu_hat <- params$estimate[2]
S_parametric <- S_FHT(sorted_times, mu, y0)
lines(sorted_times, S_parametric, col='red')



rho_mu <- function(y0, mu, sigma2, t, delta) {
  # t is vector of times, delta is vector of observed = 1, or not = 0
  log_f <- - (y0 + mu*t)/sigma2
  log_S1 <- (t/sqrt(sigma2*t)) * dnorm((mu*t + y0)/(sqrt(sigma2*t)))
  log_S2 <- (2*y0/sigma2) * exp(-2*y0*mu/sigma2) * pnorm((mu*t-y0)/sqrt(sigma2*t))
  log_S3 <- -t/sqrt(sigma2*t) * exp(-2*y0*mu/sigma2)*dnorm((mu*t - y0)/sqrt(sigma2*t))
  log_S <- log_S1 + log_S2 + log_S3
  total <- delta*log_f + (1-delta)*log_S
  return(-sum(total))
}

y0s <- seq(0.1, 3, by=0.01)
rhos <- sapply(y0s, function(y0) { rho_mu(y0, mu_hat, 1, sorted_times, sorted_observed) }, simplify='array')
plot(y0s, rhos, typ='l')