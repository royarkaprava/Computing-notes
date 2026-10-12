x = rnorm(1000, 0, 1)
p = 1 / (1 + exp(3 + 2*x))
y = rbinom(1000, size = 1, prob = p)

X <- cbind(rep(1, 1000), x)

gradientdescentfitR(X, y, 10000, 1, 0.1)