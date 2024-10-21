newtonfitR <- function(X, y, Total_itr, alpha){
  n = length(y);
  p = ncol(X);
  
  beta = rep(0, p);
  
  beta0 = beta;
  
  lossold = sum((y-X%*%beta)^2)/n;
  lossnew = lossold;
  
  y_hat = X%*%beta;
  derivbeta = rep(0, p);
  
  itr = 0;
  
  changeL = 0;
  L1 = 1;
  
  Xt = t(X);
  
  hessinv = solve(2*Xt%*%X/n); #The Hessian is 2*Xt*X/n
  
  for(itr in 1:Total_itr){
    L1 = 1;
    beta0   = beta;
    lossold = lossnew;
    y_hat = X%*%beta;
    deriv = - 2*Xt%*%(y-y_hat)/n;
    
    step = hessinv%*%deriv;
    beta = beta0 - L1*step; #Update beta
    
    lossnew = sum((y-X%*%beta)^2)/n;
    
    changeL = 0;
    
    # Implementation of damped
    while(lossnew >= lossold -alpha * L1*sum(step*deriv) ){
      changeL = changeL + 1;
      L1 = L1/2;
      beta = beta0 - L1*step; #Update beta
      
      lossnew = sum((y-X%*%beta)^2)/n;
      
      if(changeL > 20){
        beta = beta0; # Could not update, so old value reassigned
        lossnew = lossold;
        break;
      }
    }
    
    #To break the for loop as well
    if(changeL > 20){
      beta = beta0; # Could not update, so old value reassigned
      lossnew = lossold;
      break;
    }
    
  }
  
  return(beta);
}


x = rnorm(1000, 0, 1)
y = 3+2*x + rnorm(1000, 0, 1)

X <- cbind(rep(1, 1000), x)

newtonfitR(X, y, 1000, 0)