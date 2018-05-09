set.seed(1)

N <- 100
d <- 5 # beta is d + 1 dimensional
p <- 3 # gamma is p + 1 dimensional
beta <- c(4, 3, 5, -1, -2, 2) # beta_0 is first element
gamma <- c(5, 1, -3, 0)

X <- cbind(rep(1, N), matrix(rnorm(d*N), ncol = d, byrow = TRUE)) # design matrix
#X_df <- data.frame(X)

Z <- cbind(rep(1, N), matrix(rnorm(p*N), ncol = p, byrow = TRUE)) # design matrix
#Z_df <- data.frame(Z)

y0 <- exp(X %*% beta)
mu <- Z %*% gamma