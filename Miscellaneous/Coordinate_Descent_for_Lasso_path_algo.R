set.seed(1234)
N = 500
p = 30
X = scale(matrix(rnorm(N*p), ncol=p))
b = c(.5,-.5, .25,-.25, .125,-.125, rep(0, 24))
y = as.numeric(scale(X %*% b + rnorm(N, sd=.5)))

# Minimize mean((y - X %*% beta)^2) + lambda * sum(abs(beta)).
lambda_max <- 2 * max(abs(crossprod(X, y))) / N
lambda_path <- exp(seq(log(lambda_max), log(0.01 * lambda_max), length.out = 100))
beta_path <- matrix(0, p, length(lambda_path))
beta <- rep(0, p)
x2 <- colSums(X^2)

# Path algorithm: retain beta between lambdas for warm starts.
for (k in seq_along(lambda_path)) {
  change <- Inf
  itr <- 0
  while (change > 1e-7 && itr < 10000) {
    betaold <- beta
    for (i in seq_len(p)) {
      residual <- y - X[, -i, drop = FALSE] %*% beta[-i]
      betastar <- sum(X[, i] * residual) / x2[i]
      cutoff <- N * lambda_path[k] / (2 * x2[i])
      beta[i] <- sign(betastar) * max(abs(betastar) - cutoff, 0)
    }
    change <- max(abs(betaold - beta))
    itr <- itr + 1
  }
  if (change > 1e-7) warning("Coordinate descent did not converge.")
  beta_path[, k] <- beta
}

matplot(log(lambda_path), t(beta_path), type = "l", lty = 1,
        xlab = expression(log(lambda)), ylab = expression(beta))

# Compare with glmnet, which uses half the squared-error loss.
fit_path <- glmnet::glmnet(X, y, alpha = 1, intercept = FALSE,
                           standardize = FALSE, lambda = lambda_path/2,
                           thresh = 1e-12)
range(beta_path - as.matrix(fit_path$beta))
