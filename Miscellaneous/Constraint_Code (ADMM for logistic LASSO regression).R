funk <- function(par){
  Xb = X %*% par
  n <- length(y)
  pvec = 1/(1+exp(Xb));
  
  if( any(Xb > 20)){
    pvec[which(Xb > 20) ] = 0 # Assigning pvec = 0 when Xb > 20 as 1/(1+exp(Xb)) will be very small
    pvec[which(Xb > 20) ] = pvec[which(Xb > 20) ]+0.0000001; #To avoid singularity log(p)
  }
  
  if( any(Xb < - 20)){
    pvec[which(Xb < - 20) ] = 1; # Assigning pvec = 1 when Xb < - 20 as 1/(1+exp(Xb)) will be close to one
    pvec[which(Xb < - 20) ] = pvec[which(Xb < - 20) ]- 0.0000001; #To avoid singularity log(1-p)
  }
  
  loss = -sum(log(pvec) * y + log(1-pvec) * (1-y));  
  
  hval <- (par-z)
  #g_hval <- out$grad
  obj = loss + 0.5 * rho * sum(hval * hval) + sum(alpha * hval)
  #G_obj = grad_loss + (rho * hval + alpha) * g_hval
  return(obj)
}

fungradk <- function(par){
  Xb = X %*% par
  n <- length(y)
  pvec = 1/(1+exp(Xb));
  
  if( any(Xb > 20)){
    pvec[which(Xb > 20) ] = 0 # Assigning pvec = 0 when Xb > 20 as 1/(1+exp(Xb)) will be very small
    pvec[which(Xb > 20) ] = pvec[which(Xb > 20) ]+0.0000001; #To avoid singularity log(p)
  }
  
  if( any(Xb < - 20)){
    pvec[which(Xb < - 20) ] = 1; # Assigning pvec = 1 when Xb < - 20 as 1/(1+exp(Xb)) will be close to one
    pvec[which(Xb < - 20) ] = pvec[which(Xb < - 20) ]- 0.0000001; #To avoid singularity log(1-p)
  }
  
  lossder <-  t(X) %*% (y-pvec)
  
  hval <- array(par-z)
  return(array(lossder) + rho * hval + alpha)
}


p <- 100
n <- 1000
beta0 <- rnorm(p, sd=10)
beta0[11:100] <- 0
X <- matrix(rnorm(n*p), n, p)
y <- rbinom(n, 1, 1/(1+exp(X%*%beta0)))


p <- ncol(X)

beta <- rep(0,p)

z <- rnorm(p)#rep(0, p)


max_iter=100; h_tol=1e-100; rho_max=1e+100
rho=1; alpha=rep(0, p); h = 1000000 ;lambda <- 0.1

for( i in 1:max_iter){
  betanew = NULL; h_new = NULL;
  while (rho < rho_max){
    betals <- optim(beta, fn=funk,gr=fungradk, method="L-BFGS-B",control = list(maxit = 1000))$par
    
    betanew <- betals
    
    znew <- betanew + alpha/rho
    znew <- sign(znew)*(abs(znew) - lambda/rho)*(abs(znew) > lambda/rho) #Soft-thresholding
    #betanew
    
    h_new = max(abs(betanew-znew))
    
    if (h_new > 0.25 * h)
      rho = rho*10
    else
      break;
  }
  beta = betanew
  z = znew
  
  h     = h_new
  alpha = alpha +  rho * (beta-z)
  
  if (h <= h_tol || rho >= rho_max)
    break;
}



sum((beta0-beta)^2)

order(abs(beta), decreasing = T)[1:10]

#glmnet is like glm() model P(y=1) as 1/(1+exp(-Xbeta)), thus compared with '+'
sum((beta0+coef(glmnet::glmnet(X, y, family = "binomial", lambda = 0.1, intercept = F))[2:(p+1)])^2)
