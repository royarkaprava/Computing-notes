#include <RcppArmadillo.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <RcppEnsmallen.h>

// [[Rcpp::depends(RcppEnsmallen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// Define a differentiable objective function by implementing both Evaluate()
// and Gradient() separately.
class LogisticRegressionFunction
{
public:
  // Construct the object with the given the design 
  // matrix and responses.
  LogisticRegressionFunction(const arma::mat& X,
                           const arma::vec& y, const double cutoff) :
  X(X), y(y), cutoff(cutoff) { }
  
  // Return the objective function for model parameters beta.
  double Evaluate(const arma::mat& beta)
  {
    colvec Xb = X*beta;
    
    colvec pvec = 1/(1+exp(Xb));
    
    if( any(Xb > cutoff)){
      pvec.elem( find(Xb > cutoff) ).zeros(); // Assigning pvec = 0 when Xb > cutoff as 1/(1+exp(Xb)) will be very small
      pvec.elem( find(Xb > cutoff) ) += 0.0000001; //To avoid singularity log(p)
    }
    
    if( any(Xb < - cutoff)){
      pvec.elem( find(Xb < - cutoff) ).ones(); // Assigning pvec = 1 when Xb < - cutoff as 1/(1+exp(Xb)) will be close to one
      pvec.elem( find(Xb < - cutoff) ) -= 0.0000001; //To avoid singularity log(1-p)
    }
    
    
    double lossold = -sum(log(pvec) % y + log(1-pvec) % (1-y))/y.n_elem;  //% is for element-wise multiplication
    
    return lossold;
  }
  
  // Compute the gradient for model parameters beta
  void Gradient(const arma::mat& beta, arma::mat& g)
  {
    colvec Xb = X*beta;
    
    colvec pvec = 1/(1+exp(Xb));
    
    if( any(Xb > cutoff)){
      pvec.elem( find(Xb > cutoff) ).zeros(); // Assigning pvec = 0 when Xb > cutoff as 1/(1+exp(Xb)) will be very small
      pvec.elem( find(Xb > cutoff) ) += 0.0000001; //To avoid singularity log(p)
    }
    
    if( any(Xb < - cutoff)){
      pvec.elem( find(Xb < - cutoff) ).ones(); // Assigning pvec = 1 when Xb < - cutoff as 1/(1+exp(Xb)) will be close to one
      pvec.elem( find(Xb < - cutoff) ) -= 0.0000001; //To avoid singularity log(1-p)
    }
    
    g =  X.t() * (y-pvec)/y.n_elem;
  }
  
private:
  // The design matrix.
  const arma::mat& X;
  // The responses to each data point.
  const arma::vec& y;
  const double& cutoff;
};


// [[Rcpp::export]]
arma::mat Logistic_reg_lbfgs(const arma::mat& X, const arma::vec& y, const arma::mat betaini, const int maxiter, const double cutoff) {
  
  // Construct the first objective function.
  LogisticRegressionFunction lrf(X, y, cutoff);
  
  // Create the L_BFGS optimizer with default parameters.
  // The ens::L_BFGS type can be replaced with any ensmallen optimizer that can
  // handle differentiable functions.
  ens::L_BFGS lbfgs;
  
  lbfgs.MaxIterations() = maxiter;
  
  // Create a starting point for our optimization randomly.
  // The model has p parameters, so the shape is p x 1.
  arma::mat beta = betaini;//(X.n_cols, 1, arma::fill::randn);
  
  // Run the optimization
  lbfgs.Optimize(lrf, beta);
  
  return beta;
}