# event <- 1 or 0
event <- 1
censored <- 1 - event
y0_i <- exp(Xt*g)
mu_i <- Xt*b
log_f_i <- log(y0) - 0.5*log(2*pi*t^3) - (y0 + mu)^2/2*t
log_S_i <- log(pnorm(-(y0+mu)/sqrt(t)) + exp(-2*y0*mu)*pnorm((mu*t-y0)/sqrt(t)))

log_likelihood <- sum(event*log_f_i + censored*log_S_i)