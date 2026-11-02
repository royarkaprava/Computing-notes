# NONLINEAR CONJUGATE GRADIENT: SIMPLE CLASS EXAMPLES
# Base R only. Click Source in RStudio: all six plots appear on one page.
# Only the optimizer is a reusable function. Each example redefines fn and gr,
# then simulates data, fits the model, compares with R, and plots directly.

########################## Nonlinear CG #####################################
# Same notation as the attached LinearCG code, except rk = -gr(xk).
# Backtracking chooses alpha; PR+ combines the new gradient and old direction.
# tol applies to the SQUARED gradient norm: tol=1e-10 means ||gr|| <= 1e-5.

NonlinearCG <- function(fn, gr, x0, tol=1e-10, max_iter=5000,
                        method=c("PR+", "SD")){
  method <- match.arg(method)
  if(!is.numeric(x0) || !length(x0) || any(!is.finite(x0)))
    stop("x0 must be a finite numeric vector.")
  if(length(tol) != 1 || !is.finite(tol) || tol <= 0)
    stop("tol must be positive and finite.")
  if(length(max_iter) != 1 || !is.finite(max_iter) || max_iter < 1 ||
     max_iter != floor(max_iter)) stop("max_iter must be a positive integer.")
  
  xk <- as.numeric(x0)
  fk <- fn(xk)
  gk <- as.numeric(gr(xk))
  if(length(fk) != 1 || !is.finite(fk) || length(gk) != length(xk) ||
     any(!is.finite(gk))) stop("Check fn and gr at the starting value.")
  rk <- -gk
  dk <- rk
  error <- sum(rk^2)
  
  num_iter <- 0L
  alphavec <- stepvec <- betavec <- numeric(max_iter)
  fvec <- gnormvec <- numeric(max_iter + 1L)
  fvec[1] <- fk
  gnormvec[1] <- sqrt(error)
  message <- "Maximum iterations reached."
  
  while(error > tol && num_iter < max_iter){
    # If the old direction no longer points downhill, restart with -gradient.
    slope <- sum(gk*dk)
    if(!is.finite(slope) || slope > -1e-3*error){
      dk <- rk
      slope <- -error
    }
    
    # Armijo backtracking: halve alpha until the objective decreases enough.
    alpha <- 1
    accepted <- FALSE
    for(j in 1:60){
      xnew <- xk + alpha*dk
      fnew <- fn(xnew)
      if(is.finite(fnew) && fnew <= fk + 1e-4*alpha*slope){
        accepted <- TRUE
        break
      }
      alpha <- alpha/2
    }
    if(!accepted){
      message <- "Line search could not find a decreasing step."
      break
    }
    
    gnew <- as.numeric(gr(xnew))
    if(length(gnew) != length(xk) || any(!is.finite(gnew))){
      message <- "Gradient is not finite at the proposed new point."
      break
    }
    rnew <- -gnew
    errornew <- sum(rnew^2)
    if(!is.finite(errornew)){
      message <- "Squared gradient norm overflowed; rescale the problem."
      break
    }
    beta <- if(method == "PR+") max(0, sum(rnew*(rnew-rk))/error) else 0
    dnew <- rnew + beta*dk
    
    num_iter <- num_iter + 1L
    alphavec[num_iter] <- alpha
    # Compute the displacement BEFORE replacing the old direction dk.
    stepvec[num_iter] <- alpha*sqrt(sum(dk^2))
    betavec[num_iter] <- beta
    fvec[num_iter + 1L] <- fnew
    gnormvec[num_iter + 1L] <- sqrt(errornew)
    
    xk <- xnew
    fk <- fnew
    gk <- gnew
    rk <- rnew
    dk <- dnew
    error <- errornew
  }
  
  converged <- error <= tol
  if(converged) message <- "Gradient tolerance met."
  steps <- seq_len(num_iter)
  iterates <- seq_len(num_iter + 1L)
  list(sol=xk, value=fk, alpha=alphavec[steps],
       step_norm=stepvec[steps], beta=betavec[steps],
       objective=fvec[iterates], gradient_norm=gnormvec[iterates],
       num_iter=num_iter,
       converged=converged, message=message)
}

########################## Plot setup ######################################
# Each row shows one fitted model and its objective over iterations.
oldpar <- par(mfrow=c(3, 2), mar=c(4.1, 4.2, 2.5, 0.8), mgp=c(2.5, 0.8, 0))
blue <- "#126B9A"
orange <- "#C26A22"
grey <- adjustcolor("grey30", alpha.f=0.4)

########################## Example 1: least squares #########################
set.seed(2026)
n <- 250
p <- 30
X <- matrix(rnorm(n*p), n, p)
beta_true <- rnorm(p)
y <- as.numeric(X %*% beta_true + rnorm(n))

fn <- function(beta) mean((as.numeric(X %*% beta)-y)^2)/2
gr <- function(beta) as.numeric(crossprod(X, X %*% beta-y)/n)

system.time(fit_ls <- NonlinearCG(fn, gr, rep(0, p), tol=1e-12))
beta_lm <- as.numeric(lm(y ~ X-1)$coefficients)
beta_solve <- as.numeric(solve(crossprod(X), crossprod(X, y)))
cat("\nLEAST SQUARES: coefficient MSE\n")
print(c(CG=mean((fit_ls$sol-beta_true)^2),
        lm=mean((beta_lm-beta_true)^2), solve=mean((beta_solve-beta_true)^2)))
cat("Largest difference from lm:", max(abs(fit_ls$sol-beta_lm)), "\n")
cat("CG iterations:", fit_ls$num_iter, "--", fit_ls$message, "\n")

yhat <- as.numeric(X %*% fit_ls$sol)
plot(y, yhat, pch=16, col=grey, xlab="Observed response", ylab="Fitted response",
     main="1. Least-squares fit")
abline(0, 1, col=blue, lwd=2)
legend("topleft", "Identity line", col=blue, lty=1, lwd=2, bty="n")
plot(0:fit_ls$num_iter, fit_ls$objective, type="l", lwd=2, col=blue,
     xlab="Iteration", ylab="Mean squared error / 2", main="1. Objective decreases")

########################## Example 2: nonlinear regression ##################
# Mean response: a*(1-exp(-b*x)). Optimize log(a),log(b) to keep a,b positive.
set.seed(27)
x <- seq(0.05, 5, length.out=100)
a_true <- 3
b_true <- 0.8
y_growth <- a_true*(1-exp(-b_true*x)) + rnorm(length(x), sd=0.2)

fn <- function(theta){
  a <- exp(theta[1])
  b <- exp(theta[2])
  mean((a*(1-exp(-b*x))-y_growth)^2)/2
}
gr <- function(theta){
  a <- exp(theta[1])
  b <- exp(theta[2])
  mu <- a*(1-exp(-b*x))
  residual <- mu-y_growth
  c(mean(residual*mu), mean(residual*a*b*x*exp(-b*x)))
}

theta0 <- log(c(a=1, b=0.15))
system.time(fit_growth <- NonlinearCG(fn, gr, theta0, tol=1e-12))
fit_nls <- nls(y_growth ~ a*(1-exp(-b*x)), start=list(a=1, b=0.15))
cat("\nNONLINEAR REGRESSION: parameter estimates\n")
print(cbind(truth=c(a=a_true, b=b_true), CG=exp(fit_growth$sol), nls=coef(fit_nls)))
cat("CG iterations:", fit_growth$num_iter, "--", fit_growth$message, "\n")

xx <- seq(0, 5, length.out=300)
a_hat <- exp(fit_growth$sol[1])
b_hat <- exp(fit_growth$sol[2])
plot(x, y_growth, pch=16, col=grey, ylim=range(0, y_growth),
     xlab="x", ylab="Response", main="2. Nonlinear growth curve")
lines(xx, exp(theta0[1])*(1-exp(-exp(theta0[2])*xx)), col=orange, lty=3, lwd=2)
lines(xx, a_true*(1-exp(-b_true*xx)), col="black", lty=2, lwd=2)
lines(xx, a_hat*(1-exp(-b_hat*xx)), col=blue, lwd=2.5)
legend("bottomright", c("Starting curve", "True mean", "Final CG fit"),
       col=c(orange, "black", blue), lty=c(3, 2, 1), lwd=2,
       bg="white", box.col="white", cex=0.8)
plot(0:fit_growth$num_iter, fit_growth$objective, type="l", lwd=2, col=blue,
     xlab="Iteration", ylab="Mean squared error / 2", main="2. Objective decreases")

########################## Example 3: logistic regression ###################
set.seed(28)
x_logistic <- seq(-3, 3, length.out=250)
X_logistic <- cbind(1, x_logistic)
beta_true_logistic <- c(intercept=-0.5, slope=1.4)
prob_true <- plogis(as.numeric(X_logistic %*% beta_true_logistic))
y_logistic <- rbinom(length(x_logistic), 1, prob_true)

fn <- function(beta){
  eta <- as.numeric(X_logistic %*% beta)
  # Stable computation of log(1+exp(eta)).
  mean(pmax(eta, 0) + log1p(exp(-abs(eta))) - y_logistic*eta)
}
gr <- function(beta){
  prob <- plogis(as.numeric(X_logistic %*% beta))
  as.numeric(crossprod(X_logistic, prob-y_logistic)/length(y_logistic))
}

system.time(fit_logistic <- NonlinearCG(fn, gr, c(0, 0), tol=1e-12))
fit_glm <- glm(y_logistic ~ x_logistic, family=binomial())
cat("\nLOGISTIC REGRESSION: coefficient estimates\n")
print(cbind(truth=beta_true_logistic, CG=fit_logistic$sol, glm=coef(fit_glm)))
cat("CG iterations:", fit_logistic$num_iter, "--", fit_logistic$message, "\n")

plot(x_logistic, 0.02+0.96*y_logistic, pch=16, col=grey, ylim=c(0, 1),
     xlab="x", ylab="Pr(Y = 1 | x)", main="3. Logistic probability curve")
abline(h=0.5, col=orange, lty=3, lwd=2)
lines(x_logistic, prob_true, col="black", lty=2, lwd=2)
lines(x_logistic, plogis(as.numeric(X_logistic %*% fit_logistic$sol)), col=blue, lwd=2.5)
legend("topleft", c("Start: probability 0.5", "True probability", "Final CG fit"),
       col=c(orange, "black", blue), lty=c(3, 2, 1), lwd=2,
       bg="white", box.col="white", cex=0.8)
plot(0:fit_logistic$num_iter, fit_logistic$objective, type="l", lwd=2, col=blue,
     xlab="Iteration", ylab="Average negative log likelihood", main="3. Objective decreases")
par(oldpar)

# Optional: plot(fit_growth$step_norm, type="l") shows the parameter-step sizes.
# To compare with steepest descent, use method="SD" in a fit call.
# Run each example's fn, gr, and fit together: fn and gr are redefined for each model.
