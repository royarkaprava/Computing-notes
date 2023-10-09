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
colvec newtonfit(mat &X, vec &y, int Total_itr, double alpha){
  int n = y.n_elem;
  int p = X.n_cols;
  
  colvec beta = zeros<vec>(p);
  
  //colvec beta = ones<vec>(G*Mb);
  colvec beta0 = beta;
  
  double lossold = sum(pow(y-X*beta, 2));
  double lossnew = lossold;
  
  colvec y_hat = X*beta;
  colvec step = zeros<vec>(p);
  colvec deriv = zeros<vec>(p);
  
  int itr = 0;
  
  int changeL = 0;
  double L1 = 1;
  
  mat Xt = X.t();
  
  mat hessinv = inv(2*Xt*X/n); //The Hessian is 2*Xt*X/n
  
  for(itr = 0; itr < Total_itr; itr++){
    L1 = 1;
    beta0   = beta;
    lossold = lossnew;
    y_hat = X*beta;
    deriv = - 2*Xt*(y-y_hat)/n;
    
    step = hessinv*deriv;
    beta = beta0 - L1*step; //Update beta
    
    lossnew = sum(pow(y-X*beta, 2));
    
    changeL = 0;
    
    // Implementation of Section 5.2 
    while(lossnew >= lossold -alpha * L1*sum(step%deriv) ){
      changeL = changeL + 1;
      L1 = L1/2;
      beta = beta0 - L1*step; //Update beta
      
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