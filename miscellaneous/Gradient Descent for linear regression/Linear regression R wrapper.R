Rcpp::sourceCpp("~/Course to teach/Computing/Gradient descent linear regression/Linear_gradient_descent.cpp")
Rcpp::sourceCpp("~/Course to teach/Computing/Gradient descent linear regression/Newton_for_linear.cpp")

x = rnorm(1000, 0, 1)
y = 3+2*x + rnorm(1000, 0, 1)

X <- cbind(rep(1, 1000), x)
betaGD = gradientdescentfit(X, y, 1000, 0.01, 0)

betaNR = newtonfit(X, y, 1000, 0)