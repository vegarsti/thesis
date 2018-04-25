set.seed(1)

# 1D example of L2Boost

N <- 1000

x <- runif(n = N, min = 0, max = 10)
eps <- rnorm(n = N, mean = 0, sd = 10)
f <- function(x) 10*sin(x) + 2*x
y <- f(x) + eps
plot(x, y)
points(x,f(x), col='red')

# Algorithm

M <- 100+1
nu <- 0.1

errors <- rep(NA, M+1)
u <- matrix(nrow = M+1, ncol = N)
F_ <- matrix(nrow = M+1, ncol = N)
F_[1, ] <- rep(mean(y), N)

F_at <- function(F_, m) {
  # Calculate function values until m
  y_hat <- rep(0, N)
  for (i in 1:N) {
    y_hat[i] <- sum(F_[1:m, i])
  }
  return(y_hat)
}

u[1, ] <- y

for (m in 2:M) {
  y_hat_m <- F_at(F_, m-1)
  u[m, ] <- (y - y_hat_m)
  
  # lm 
  h <- lm(u[m, ] ~ x)
  beta_0 <- h$coefficients[1]
  beta_1 <- h$coefficients[2]
  F_[m, ] <- beta_0 + nu*beta_1*x
  
  # smoothing spline
  #s <- smooth.spline(u[m, ]~x, df=2.5)
  #F_[m, ] <- predict(s, x)$y
}
y_hat_M <- F_at(F_, M)
points(x, y_hat_M, col='blue')

a <- smooth.spline(y~x, df=10)
y_hat <- predict(a, x)$y
points(x, y_hat, col='yellow')