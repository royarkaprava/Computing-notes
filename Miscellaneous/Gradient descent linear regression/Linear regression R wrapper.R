Rcpp::sourceCpp("~/Course to teach/Computing/Gradient descent linear regression/Linear_gradient_descent.cpp")

x = rnorm(1000, 0, 1)
y = 3+2*x + rnorm(1000, 0, 1)

X <- cbind(rep(1, 1000), x)

gradientdescentfitR(X, y, 1000, 0.01, 0)

##Using the Rcpp code
betaGD = gradientdescentfit(X, y, 1000, 0.01, 0)