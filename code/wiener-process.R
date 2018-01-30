w <- rep(NA, 1000)
w[1] <- 0

for (i in 2:1000) {
  w[i] <- w[i-1] + rnorm(n=1, mean=0, sd=1)
}

plot.ts(w)