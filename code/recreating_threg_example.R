library(threg)
data("lkr")
lkr$f.treatment2=factor(lkr$treatment2)
fit <- threg(Surv(weeks, relapse) ~ f.treatment2|f.treatment2, data=lkr)

library(dplyr)
#tbl <- select(.data=lkr, weeks, relapse, f.treatment2)
tbl <- data.frame(lkr$weeks, lkr$relapse, lkr$f.treatment2,
                  row_names = c("weeks", "relapse", "f.treatment2"))
names(tbl)[names(tbl) == "lkr.weeks"] <- "weeks"
names(tbl)[names(tbl) == "lkr.relapse"] <- "relapse"
names(tbl)[names(tbl) == "lkr.f.treatment2"] <- "f.treatment2"

n <- dim(tbl)[1]

to_optimize <- function(params) {
  total_loglikelihood <- 0
  
  gamma <- params[1:2]
  beta <- params[3:4]
  
  for (i in 1:n) {
    tbl_i <- tbl[i, ]
    event <- tbl_i$relapse
    t_i <- tbl_i$weeks
    is_treated <- as.integer(tbl_i$f.treatment2)-1
    X_i <- c(1, is_treated)
    y0_i <- exp(sum(gamma*X_i))
    mu_i <- sum(beta*X_i)
    log_f_i <- log(y0_i) - 0.5*log(2*pi*t_i^3) - ((y0_i + mu_i*t_i)^2)/(2*t_i)
    log_S_i <- log(1 - pnorm(-(y0_i+mu_i*t_i)/sqrt(t_i)) - exp(-2*y0_i*mu_i)*pnorm((mu_i*t_i-y0_i)/sqrt(t_i)))
    loglik_i <- event*log_f_i + (1 - event)*log_S_i
    total_loglikelihood <- total_loglikelihood + loglik_i
  }
  return(-total_loglikelihood)
}

params_from_threg <- c(2.0098, -1.2739, -0.5886, 0.5886)
threg_value <- -to_optimize(params_from_threg)

initial_params <- c(1, 1, 1, 1)

best <- nlm(to_optimize, initial_params)
params_from_best <- best$estimate
best_value <- -best$minimum
print(best_value)
print(threg_value)