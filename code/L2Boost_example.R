set.seed(1)

N <- 100
x1 <- runif(n = N, min = 0, max = 10)
x2 <- runif(n = N, min = 0, max = 10)

eps <- rnorm(n = N, mean = 0, sd = 2)
y <- sin(x1) + 0.5*sin(x2) + eps

# L2 boost

M <- 100 + 1
f <- matrix(nrow = M, ncol = 3)
f[1, 1] <- mean(y)
f[1, 2] <- 0
f[1, 3] <- 0

nu <- 0.1

for (m in 2:M) {
  y_hat_m <- rep(0, N)
  for (k in 1:(m-1)) {
    y_hat_m <- y_hat_m + f[k, 2]*x1 + f[k, 3]*x2
  }
  u <- y - y_hat_m
  g <- lm(u ~ x1 + x2)
  f[m, ] <- f[m-1, ]
  f[m, 2] <- f[m, 2] + nu*g$coefficients[2]
  f[m, 3] <- f[m, 3] + nu*g$coefficients[3]
}