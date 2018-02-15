library(threg)
data("lkr")
lkr$f.treatment2=factor(lkr$treatment2)
# head(lkr)
fit <- threg(Surv(weeks, relapse) ~ f.treatment2|f.treatment2, data=lkr)
fit