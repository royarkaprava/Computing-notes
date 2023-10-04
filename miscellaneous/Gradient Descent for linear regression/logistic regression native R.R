logisticGDfitR <- function(X, y, Total_itr, L, alpha){
  n = length(y);
  p = ncol(X);
  
  beta = rep(0, p);
  
  beta0 = beta;
  
  Xt = t(X);
  Xb = X%*%beta;
  
  pvec = 1/(1+exp(Xb));
  
  if( any(Xb > 200)){
    pvec[which(Xb > 200) ] = 0 # Assigning pvec = 0 when Xb > 200 as 1/(1+exp(Xb)) will be very small
    pvec[which(Xb > 200) ] = pvec[which(Xb > 200) ]+0.0000001; #To avoid singularity log(p)
  }
  
  if( any(Xb < - 200)){
    pvec[which(Xb < - 200) ] = 1; # Assigning pvec = 1 when Xb < - 200 as 1/(1+exp(Xb)) will be close to one
    pvec[which(Xb < - 200) ] = pvec[which(Xb < - 200) ]- 0.0000001; #To avoid singularity log(1-p)
  }
  
  
  lossold = -sum(log(pvec) * y + log(1-pvec) * (1-y))/n;  
  lossnew = lossold;
  
  
  derivbeta = rep(0, p);
  
  itr = 0;
  loss = rep(0, Total_itr);
  
  changeL = 0;
  L1 = L;
  
  losslist =  rep(0, Total_itr);
  
  for(itr in 1:Total_itr){
    if(changeL > 20){
      beta = beta0; # Could not update, so old value reassigned
      lossnew = lossold;
      break;
    }
    lossold = lossnew;
    L1 = L;
    beta0 = beta;
    
    derivbeta =  Xt %*% (y-pvec)/n;  
    beta = beta0 - L1*derivbeta; #Update beta
    Xb = X%*%beta;
    
    pvec = 1/(1+exp(Xb)); 
    
    if( any(Xb > 200)){
      pvec[which(Xb > 200) ] =0 # Assigning pvec = 0 when Xb > 200 as 1/(1+exp(Xb)) will be very small
      pvec[which(Xb > 200) ] = pvec[which(Xb > 200) ]+0.0000001; #To avoid singularity log(p)
    }
    
    if( any(Xb < - 200)){
      pvec[which(Xb < - 200) ] = 1; # Assigning pvec = 1 when Xb < - 200 as 1/(1+exp(Xb)) will be close to one
      pvec[which(Xb < - 200) ] = pvec[which(Xb < - 200) ]- 0.0000001; #To avoid singularity log(1-p)
    }
    
    lossnew = -sum(log(pvec) * y + log(1-pvec) * (1-y))/n; 
    
    changeL = 0;
    
    
    while(lossnew >= lossold - alpha * L1*sum(derivbeta*derivbeta) ){
      changeL = changeL + 1;
      L1 = L1/2;
      beta = beta0 - L1*derivbeta; #Update beta
      Xb = X%*%beta;
      
      pvec = 1/(1+exp(Xb)); 
      
      if( any(Xb > 200)){
        pvec[which(Xb > 200) ] =0 # Assigning pvec = 0 when Xb > 200 as 1/(1+exp(Xb)) will be very small
        pvec[which(Xb > 200) ] = pvec[which(Xb > 200) ]+0.0000001; #To avoid singularity log(p)
      }
      
      if( any(Xb < - 200)){
        pvec[which(Xb < - 200) ] = 1; # Assigning pvec = 1 when Xb < - 200 as 1/(1+exp(Xb)) will be close to one
        pvec[which(Xb < - 200) ] = pvec[which(Xb < - 200) ]- 0.0000001; #To avoid singularity log(1-p)
      }
      
      lossnew = -sum(log(pvec) * y + log(1-pvec) * (1-y))/n; 
      
      if(changeL > 20){
        beta = beta0; # Could not update, so old value reassigned
        lossnew = lossold;
        break;
      }
    }
    
    if(changeL > 20){
      beta = beta0; # Could not update, so old value reassigned
      lossnew = lossold;
      break;
    }
    
    losslist[itr] = lossnew; 
  }
  
  
  ret = list(beta=beta, loss=losslist);
  
  return(ret)
}
