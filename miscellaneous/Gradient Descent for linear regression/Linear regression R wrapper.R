Rcpp::sourceCpp("Linear_gradient_descent.cpp")
source("linear regression native R.R")

x = rnorm(1000, 0, 1)
y = 3+2*x + rnorm(1000, 0, 1)

X <- cbind(rep(1, 1000), x)
system.time(betaGD <- gradientdescentfit(X, y, 1000, 0.01, 0))

system.time(betaGDR <- gradientdescentfitR(X, y, 1000, 0.01, 0))
