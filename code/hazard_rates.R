f <- function(t, c, mu, sigma) {
  return(c/(sigma*sqrt(2*pi))*t^(-3/2)*exp(-(c-mu*t)^2/(2*sigma^2*t)))
}

ts <- seq(from=0, to=4, by=0.01)
y <- f(ts, c=0.2, mu=1, sigma=1)
plot.ts(y, ylim=c(0, 2))
for (c in c(1, 3)) {
  y <- f(ts, c=c, mu=1, sigma=1)
  lines(y)
}