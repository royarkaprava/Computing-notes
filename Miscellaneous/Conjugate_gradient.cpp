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
colvec LinearCGc(mat &A, colvec &b, colvec &x0, double tol=1e-5){
  colvec xk = x0;
  colvec rk = b - A * xk;
  colvec dk = rk;
  double error = sum(rk%rk);
  
  colvec adk = A * dk;
  double rk2old = error;
  
  double alpha = rk2old / sum(dk%adk);
  
  double beta = sum(rk%rk) / rk2old;
  
  int n = A.n_rows;
  
  int num_iter = 0;
  
  while(error > tol){
    adk = A * dk;
    rk2old = error;
    
    alpha = rk2old / sum(dk%adk);
    xk = xk + alpha * dk;
    rk = rk - alpha * adk;
    beta = sum(rk%rk) / rk2old;
    dk = rk + beta * dk;
    
    num_iter = num_iter + 1;
    
    if(num_iter > n){break;}
    error = sum(rk%rk);
  }

  return(xk);
}
