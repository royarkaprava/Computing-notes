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
colvec PreLinearCGc(colvec &minv, mat &A, colvec &b, colvec &x0, double tol=1e-5){
  colvec xk = x0;
  colvec rk = b - A * xk;
  colvec zk = minv % rk;
  colvec dk = zk;
  double error = sum(rk%rk);
  
  colvec adk = A * dk;
  double rkzk = sum(rk%zk);
  
  double alpha = rkzk / sum(dk%adk);
  
  double beta = sum(rk%rk) / rkzk;
  
  int n = A.n_rows;
  
  int num_iter = 0;
  
  while(error > tol){
    adk = A * dk;
    rkzk = sum(rk%zk);
    
    alpha = rkzk / sum(dk%adk);
    xk = xk + alpha * dk;
    rk = rk - alpha * adk;
    zk = minv % rk;
    
    beta = sum(rk%zk) / rkzk;
    dk = zk + beta * dk;
    
    num_iter = num_iter + 1;
    
    if(num_iter > n){break;}
    error = sum(rk%rk);
  }
  
  return(xk);
}