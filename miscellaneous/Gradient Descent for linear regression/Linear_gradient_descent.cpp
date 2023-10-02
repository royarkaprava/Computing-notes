#include <RcppArmadillo.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include<omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;


// The part below is only needed if called in R.

// [[Rcpp::export]]
colvec gradientdescentfit(mat &X, vec &y, int Total_itr, double L, double alpha){
  int n = y.n_elem;
  int p = X.n_cols;
  
  colvec beta = zeros<vec>(p);
  
  //colvec beta = ones<vec>(G*Mb);
  colvec beta0 = beta;
  
  double lossold = sum(pow(y-X*beta, 2));
  double lossnew = lossold;
  
  colvec y_hat = X*beta;
  colvec derivbeta = zeros<vec>(p);
 
  int itr = 0;
  
  int changeL = 0;
  double L1 = L;
  
  mat Xt = X.t();
  
  for(itr = 0; itr < Total_itr; itr++){
    L1 = L;
    beta0   = beta;
    lossold = lossnew;
    y_hat = X*beta;
    
    derivbeta = - 2*Xt*(y-y_hat)/n;
    beta = beta0 - L1*derivbeta; //Update beta
    
    lossnew = sum(pow(y-X*beta, 2));
    
    changeL = 0;
    
    // Implementation of Section 5.2 
    while(lossnew >= lossold -alpha * L1*sum(derivbeta%derivbeta) ){
      changeL = changeL + 1;
      L1 = L1/2;
      beta = beta0 - L1*derivbeta; //Update beta
      
      lossnew = sum(pow(y-X*beta, 2));
      
      if(changeL > 20){
        beta = beta0; // Could not update, so old value reassigned
        lossnew = lossold;
        break;
      }
    }
    
    //To break the for loop as well
    if(changeL > 20){
      beta = beta0; // Could not update, so old value reassigned
      lossnew = lossold;
      break;
    }
    
  }
  
  return beta;
}