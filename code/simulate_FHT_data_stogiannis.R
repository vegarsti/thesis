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
    t <- censored_survival_times
    log_f <- log(y0) - 0.5*log(2*pi*t^3) - ((y0 + mu*t)^2)/(2*t)
    log_S <- log(pnorm((y0+mu*t)/sqrt(t)) - exp(-2*y0*mu)*pnorm((mu*t-y0)/sqrt(t)))
    loglik <- sum(observed*log_f) + sum((1 - observed)*log_S)
    total_loglikelihood <- sum(loglik)
    return(-total_loglikelihood)
  }
  
  N_params <- p+d+2
  initial_params <- runif(N_params, min=1, max=2)
  fit <- nlm(loglikfunc, initial_params)
  fit
}
### SIMULATION EXPERIMENT
set.seed(1)
# Same scheme as in "Comparing first hitting time and proportional hazards regression models" (Stogiannis, 2010)
N <- 100
d <- 1
p <- 1
x1 <- c(rep(1, 50), rep(0, 50))
x2 <- rbinom(N, size=1, p=0.5)
x3 <- runif(N, min=0, max=1)
x4 <- rnorm(2, 1)
X <- cbind(x1, x2, x3, x4)
beta <- c(0, 2, -1, 2, -1)
gamma <- c(0, 0.1, 0.2, 0.3, 0.4)
#Z <- cbind(c(rep(1, 200), rep(4, 300), rep(10, 500)))
y_intercepts <- rep(1, N)
X_design_matrix <- cbind(y_intercepts, X)
Z_design_matrix <- cbind(y_intercepts, X)
y0 <- exp(X_design_matrix %*% beta)
mu <- Z_design_matrix %*% gamma
sigma_2 <- 1
mu_IG <- y0/abs(mu)
lambda_IG <- (y0/sigma_2)^2
survival_times <- rinvgauss(N, mu_IG, lambda_IG)
W <- 0.7 # censoring time
censored_survival_times <- survival_times
censored_survival_times[is.na(survival_times)] <- W # censor!
observed <- ifelse(survival_times < W, 1, 0)
times <- pmin(survival_times, W)


fit <- run_nlm(X, X, censored_survival_times, observed)
d <- length(beta)
p <- length(gamma)
beta_hat <- fit$estimate[1:1+d]
gamma_hat <- fit$estimate[d+2:2+d+p]
beta
beta_hat
gamma
gamma_hat