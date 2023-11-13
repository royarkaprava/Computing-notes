######################data generation###################################################
p <- 10
n <- 1000
beta0 <- rnorm(p, sd=10)
Xf <- matrix(rnorm(n*10), n, p)
yf <- rbinom(n, 1, 1/(1+exp(Xf%*%beta0)))


cutoff = 20

#################################################################################################
funf <- function(par){
  Xb = Xf %*% par
  n <- length(yf)
  pvec = 1/(1+exp(Xb));
  
  if( any(Xb > cutoff)){
    pvec[which(Xb > cutoff) ] = 0 # Assigning pvec = 0 when Xb > cutoff as 1/(1+exp(Xb)) will be very small
    pvec[which(Xb > cutoff) ] = pvec[which(Xb > cutoff) ]+0.0000001; #To avoid singularity log(p)
  }
  
  if( any(Xb < - cutoff)){
    pvec[which(Xb < - cutoff) ] = 1; # Assigning pvec = 1 when Xb < - cutoff0 as 1/(1+exp(Xb)) will be close to one
    pvec[which(Xb < - cutoff) ] = pvec[which(Xb < - cutoff) ]- 0.0000001; #To avoid singularity log(1-p)
  }
  
  loss = -sum(log(pvec) * yf + log(1-pvec) * (1-yf))/n;  
  
  obj = loss + sum(theta * abs(par))
  #G_obj = grad_loss + (rho * hval + alpha) * g_hval
  return(obj)
}

fungradf <- function(par){
  Xb = Xf %*% par
  n <- length(yf)
  pvec = 1/(1+exp(Xb));
  
  if( any(Xb > cutoff)){
    pvec[which(Xb > cutoff) ] = 0 # Assigning pvec = 0 when Xb > cutoff as 1/(1+exp(Xb)) will be very small
    pvec[which(Xb > cutoff) ] = pvec[which(Xb > cutoff) ]+0.0000001; #To avoid singularity log(p)
  }
  
  if( any(Xb < - cutoff)){
    pvec[which(Xb < - cutoff) ] = 1; # Assigning pvec = 1 when Xb < - cutoff as 1/(1+exp(Xb)) will be close to one
    pvec[which(Xb < - cutoff) ] = pvec[which(Xb < - cutoff) ]- 0.0000001; #To avoid singularity log(1-p)
  }
  
  lossder <-  t(Xf) %*% (yf-pvec)/n
  
  return(lossder +theta*sign(par))
}
betaf <- rep(0, p)

maxiter <- 10000

############Using R's optim###############################

betaf <- optim(betaf, fn=funf,gr=fungradf, method="BFGS",control = list(maxit = maxiter))$par

sum((beta0-betaf)^2)

betaf <- optim(betaf, fn=funf,gr=fungradf, method="L-BFGS-B",control = list(maxit = maxiter))$par

sum((beta0-betaf)^2)

betaf <- optim(betaf, fn=funf,gr=fungradf, method="L-BFGS-B",
               lower = rep(-100,p), upper = rep(100,p),control = list(maxit = maxiter, trace=3))$par

sum((beta0-betaf)^2)

############Using RcppEnsmallen's optimization tool###############################
Rcpp::sourceCpp("Logistic_regression_Ensmallen.cpp")

coefs_lbfgs <- Logistic_reg_lbfgs(Xf, yf, matrix(betaf, p, 1), maxiter, cutoff)

sum((beta0-coefs_lbfgs)^2)

#################Using glm's exact newton solution##########
#glmnet is like glm() model P(y=1) as 1/(1+exp(-Xbeta)), thus compared with '+'
sum((beta0+coef(glm(yf~Xf-1, family = "binomial")))^2)
