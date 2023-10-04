#include <RcppArmadillo.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <omp.h>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// [[Rcpp::export]]
Rcpp::List logisticGDfit(mat &X, vec &y, int Total_itr, double L, double alpha){
  int n = y.n_elem;
  int p = X.n_cols;
  
  colvec beta = zeros<vec>(p);
  
  //colvec beta = ones<vec>(G*Mb);
  colvec beta0 = beta;
  
  mat Xt = X.t();
  colvec Xb = X*beta;
  
  colvec pvec = 1/(1+exp(Xb));
  
  if( any(Xb > 20)){
    pvec.elem( find(Xb > 20) ).zeros(); // Assigning pvec = 0 when Xb > 200 as 1/(1+exp(Xb)) will be very small
    pvec.elem( find(Xb > 20) ) += 0.0000001; //To avoid singularity log(p)
  }
  
  if( any(Xb < - 20)){
    pvec.elem( find(Xb < - 20) ).ones(); // Assigning pvec = 1 when Xb < - 200 as 1/(1+exp(Xb)) will be close to one
    pvec.elem( find(Xb < - 20) ) -= 0.0000001; //To avoid singularity log(1-p)
  }
  
  
  double lossold = -sum(log(pvec) % y + log(1-pvec) % (1-y))/n;  //% is for element-wise multiplication
  double lossnew = lossold;
  
  //colvec y_hat = X*beta;
  colvec derivbeta = zeros<vec>(p);
  
  int itr = 0;
  colvec loss = zeros<vec>(Total_itr);
  
  int changeL = 0;
  double L1 = L;
  
  colvec losslist =  zeros<vec>(Total_itr);
  
  for(itr = 0; itr < Total_itr; itr++){
    if(changeL > 20){
      beta = beta0; // Could not update, so old value reassigned
      lossnew = lossold;
      break;
    }
    lossold = lossnew;
    L1 = L;
    beta0 = beta;
    
    derivbeta =  Xt * (y-pvec)/n;  
    beta = beta0 - L1*derivbeta; //Update beta
    Xb = X*beta;
    
    pvec = 1/(1+exp(Xb)); 
    
    if( any(Xb > 20)){
      pvec.elem( find(Xb > 20) ).zeros(); // Assigning pvec = 0 when Xb > 20 as 1/(1+exp(Xb)) will be very small
      pvec.elem( find(Xb > 20) ) += 0.0000001; //To avoid singularity log(p)
    }
    
    if( any(Xb < - 20)){
      pvec.elem( find(Xb < - 20) ).ones(); // Assigning pvec = 1 when Xb < - 20 as 1/(1+exp(Xb)) will be close to one
      pvec.elem( find(Xb < - 20) ) -= 0.0000001; //To avoid singularity log(1-p)
    }
    
    lossnew = -sum(log(pvec) % y + log(1-pvec) % (1-y))/n; //% is for element-wise multiplication
    
    changeL = 0;
    
    // Implementation of Section 5.2 
    while(lossnew >= lossold - alpha * L1*sum(derivbeta%derivbeta) ){
      changeL = changeL + 1;
      L1 = L1/2;
      beta = beta0 - L1*derivbeta; //Update beta
      Xb = X*beta;
        
      pvec = 1/(1+exp(Xb)); 
      
      if( any(Xb > 20)){
        pvec.elem( find(Xb > 20) ).zeros(); // Assigning pvec = 0 when Xb > 20 as 1/(1+exp(Xb)) will be very small
        pvec.elem( find(Xb > 20) ) += 0.0000001; //To avoid singularity in log(p)
      }
      
      if( any(Xb < - 20)){
        pvec.elem( find(Xb < - 20) ).ones(); // Assigning pvec = 1 when Xb < - 20 as 1/(1+exp(Xb)) will be close to one
        pvec.elem( find(Xb < - 20) ) -= 0.0000001; //To avoid singularity in log(1-p)
      }
      
      lossnew = -sum(log(pvec) % y + log(1-pvec) % (1-y))/n;
      
      if(changeL > 20){
        beta = beta0; // Could not update, so old value reassigned
        lossnew = lossold;
        break;
      }
    }
    
    if(changeL > 20){
      beta = beta0; // Could not update, so old value reassigned
      lossnew = lossold;
      break;
    }
     
    losslist(itr) = lossnew; 
  }
  
  Rcpp::List ret;
  ret["beta"] = beta;
  ret["loss"] = losslist;
  
  return ret;
}
