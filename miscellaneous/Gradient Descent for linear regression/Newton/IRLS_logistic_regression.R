############Code Courtesy: Weijia Jin 

##This code looks much better than what I was planning to write

getp = function(X,beta){
  p = 1/(1+exp(X%*%beta))
  return (p[,1])
}
getL = function(Y,X,beta){
  n = length(Y)
  phat = getp(X,beta)
  logL = sum(Y*log(phat)+(1-Y)*log(1-phat))
  -logL/n
}
getderive = function(Y,X,beta){
  derive = 0
  n = length(Y)
  phat = getp(X,beta)
  derive = colSums(X*(Y-phat))
  derive
}


IRLSfitR <- function(X, y, Total_itr){
  n = length(y);
  p = ncol(X);
  beta = rep(0, p);
  beta0 = beta;
  lossold = getL(y,X,beta0);
  lossnew = lossold;
  itr = 0;
  Xt= t(X)
  lossvec = c()
  
  for(itr in 1:Total_itr){
    lossvec[itr] = lossnew
    lossold = lossnew
    beta0 = beta
    phat = getp(X,beta0)
    W = diag(phat)
    Winv = diag(1/phat)
    
    z = Winv%*%(y-phat)
    u = X%*%beta0-z
    beta = solve(Xt%*%W%*%X)%*%Xt%*%W%*%u
    lossnew = getL(y,X,beta);
    #if (mean( (beta-beta0)^2/beta0^2)<0.0000001){break} # if convergence then break
  }

  return(list(beta=beta, loss=lossvec))
}

x = rnorm(1000, 0, 1)
y = rbinom(1000, 1, 1/(1+exp(-3+2*x)))

X <- cbind(rep(1, 1000), x)

out = IRLSfitR(X, y, 10000)
out$beta
plot(out$loss[which(out$loss>0)], type='l')

#Compare with glm
summary(glm(y~x, family = "binomial"))$coefficients  #The sign is opposite as R sets  pi = 1/(1+exp(-xbeta))
