t <- seq(from=0, to=20, by=0.01)

f <- function(t, y0, mu, sigma2) {
  return((y0/sqrt(2*pi*sigma2*t^3))*exp(-((y0+mu*t)^2/(2*sigma2*t))))
}

sigma2 <- 1
mu <- -1
y1 <- f(t, y0=2, mu=mu, sigma2=sigma2)
y2 <- f(t, y0=5, mu=mu, sigma2=sigma2)
y3 <- f(t, y0=10, mu=mu, sigma2=sigma2)

colors <- c("tomato2", "whitesmoke", "wheat1")

plot(t, y1, type="l", col=colors[1])
lines(t, y2, col=colors[2])
lines(t, y3, col=colors[3])
legend('topright', title="y0", legend=c(2, 5, 10), col=colors, lty=1, cex=0.8)