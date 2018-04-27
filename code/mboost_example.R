library("mboost")
data("bodyfat", package="TH.data")

# Reproduce from paper
lm1 <- lm(DEXfat ~ hipcirc + kneebreadth + anthro3a, data=bodyfat)
coef(lm1)

# Estimate same by glmboost
glm1 <- glmboost(DEXfat ~ hipcirc + kneebreadth + anthro3a, data = bodyfat)
coef(glm1, off2int=TRUE)

glm2 <- glmboost(DEXfat ~ ., data = bodyfat)
coef(glm2)

preds <- names(bodyfat[, names(bodyfat) != "DEXfat"])
plot(glm2, ylim=range(coef(glm2, which=preds)))
