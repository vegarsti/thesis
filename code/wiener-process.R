#set.seed(1)

N <- 200
t <- 1:N

M <- 10 # no. of individuals
y0 <- rep(10, M) #runif(n = M, min=10, max=20)
mu <- -0.1

y <- matrix(nrow = N, ncol = M)
delta <- rep(0, M)
y[1, ] <- y0

for (i in 2:N) {
  for (j in 1:M) {
    y[i, j] <- y[i-1, j] + rnorm(n = 1, mean = mu, sd = 1)
    if (y[i-1, j] > 0) {
      if (y[i, j] < 0) {
        delta[j] <- i
      }
    }
  }
}

plot(t, y[, 1], type='l', ylim = c(-50, 50))
abline(a = 0, b = 0)
lines(t, y[, 2], col = 'red')
lines(t, y[, 3], col = 'blue')
lines(t, y[, 4], col = 'gray')
lines(t, y[, 5], col = 'darkorange')
lines(t, y[, 6], col = 'brown')
lines(t, y[, 7], col = 'brown')
lines(t, y[, 8], col = 'darkorchid')
print(delta)
print(mean(delta))