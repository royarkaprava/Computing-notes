gradientdescentfitR <- function(X, y, Total_itr, L, alpha){
  n = length(y);
  p = ncol(X);
  
  beta = rep(0, p);
  
  beta0 = beta;
  
  lossold = sum((y-X%*%beta)^2);
  lossnew = lossold;
  
  y_hat = X%*%beta;
  derivbeta = rep(0, p);
  
  itr = 0;
  
  changeL = 0;
  L1 = L;
  
  Xt = t(X);
  
  for(itr in 1:Total_itr){
    L1 = L;
    beta0   = beta;
    lossold = lossnew;
    y_hat = X%*%beta;
    
    derivbeta = - 2*Xt%*%(y-y_hat)/n;
    beta = beta0 - L1*derivbeta; 
    
    lossnew = sum((y-X%*%beta)^2);
    
    changeL = 0;
    
    # Implementation of Section 5.2 
    while(lossnew >= lossold -alpha * L1*sum(derivbeta^2) ){
      changeL = changeL + 1;
      L1 = L1/2;
      beta = beta0 - L1*derivbeta; #Update beta
      
      lossnew = sum((y-X%*%beta)^2);
      
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