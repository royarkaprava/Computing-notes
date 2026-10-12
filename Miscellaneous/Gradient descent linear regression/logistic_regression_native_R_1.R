# Logistic regression by gradient descent with Armijo backtracking.
# PDF convention: p = 1 / (1 + exp(X %*% beta)).
# X must include a column of ones if an intercept is wanted; y must be 0/1.
# Set link_sign = 1 for the usual glm() convention, p = plogis(X %*% beta).
# The original function name and first five arguments are retained.
gradientdescentfitR <- function(X, y, Total_itr, L, alpha,
                               link_sign = -1, tol = 1e-8){
  X = as.matrix(X);
  y = as.numeric(y);
  n = length(y);
  p = ncol(X);

  stopifnot(is.numeric(X), n > 0, p > 0, nrow(X) == n,
            all(is.finite(X)), all(is.finite(y)), all(y %in% c(0, 1)),
            length(Total_itr) == 1, is.finite(Total_itr),
            Total_itr >= 1, Total_itr == floor(Total_itr),
            length(L) == 1, is.finite(L), L > 0,
            length(alpha) == 1, is.finite(alpha), alpha > 0, alpha < 1,
            length(link_sign) == 1, link_sign %in% c(-1, 1),
            length(tol) == 1, is.finite(tol), tol > 0);

  beta = rep(0, p);
  Xt = t(X);

  # Average negative log-likelihood, evaluated without log(0) or exp overflow.
  # eta is the log odds; for the PDF eta = -X %*% beta.
  loss = function(beta){
    eta = link_sign * drop(X %*% beta);
    if(any(!is.finite(eta))) return(Inf);
    mean(ifelse(y == 1, pmax(-eta, 0), pmax(eta, 0)) +
           log1p(exp(-abs(eta))));
  };

  lossnew = loss(beta);
  converged = FALSE;
  line_search_failed = FALSE;

  for(itr in seq_len(Total_itr)){
    L1 = L;
    beta0 = beta;
    lossold = lossnew;
    p_hat = plogis(link_sign * drop(X %*% beta));

    # PDF convention: gradient = X^T (y - p_hat) / n.
    derivbeta = link_sign * drop(Xt %*% (p_hat - y)) / n;
    if(max(abs(derivbeta)) <= tol){
      converged = TRUE;
      break;
    }

    beta = beta0 - L1 * derivbeta;
    lossnew = loss(beta);
    changeL = 0;

    # Use the same average loss as the gradient in the Armijo condition.
    while(lossnew > lossold - alpha * L1 * sum(derivbeta^2)){
      changeL = changeL + 1;
      if(changeL > 50){
        beta = beta0;
        lossnew = lossold;
        line_search_failed = TRUE;
        break;
      }
      L1 = L1 / 2;
      beta = beta0 - L1 * derivbeta;
      lossnew = loss(beta);
    }

    if(line_search_failed) break;
  }

  # Check the final update as well.
  p_hat = plogis(link_sign * drop(X %*% beta));
  derivbeta = link_sign * drop(Xt %*% (p_hat - y)) / n;
  converged = max(abs(derivbeta)) <= tol;

  if(line_search_failed){
    warning("Backtracking failed; returning the last accepted coefficients.");
  } else if(!converged){
    warning("Iteration limit reached before convergence; increase Total_itr. ",
            "Under separation a finite maximum-likelihood estimate may not exist.");
  }

  return(beta);
}

# Example (uncomment to run):
# set.seed(123)
# X = cbind(Intercept = 1, x1 = rnorm(500), x2 = rnorm(500));
# beta_true = c(0.5, -1, 0.8);
# y = rbinom(nrow(X), 1, plogis(-drop(X %*% beta_true)));
# beta_hat = gradientdescentfitR(X, y, Total_itr = 10000, L = 1, alpha = 0.1);
# p_hat = plogis(-drop(X %*% beta_hat));
# glm_beta = coef(glm(y ~ X - 1, family = binomial()));
# cbind(PDF_convention = beta_hat, equivalent_glm = -glm_beta);

