set.seed(1234)
N = 500
p = 30
X = scale(matrix(rnorm(N*p), ncol=p))
b = c(.5,-.5, .25,-.25, .125,-.125, rep(0, 24))
y = scale(X %*% b + rnorm(N, sd=.5))

#Lambda's are sorted in decreasing order

lambda =0.1

#Since LASSO is applied for high dimensional data, the common initialization technique
#is usually to set all to 0 for reducing computational cost in initialization.

beta <- betastar <- rep(0, p)


  for(itr in 1:10000){
    #Old beta should be stored outside of the updating loop for beta
    beta0old = beta
    for (i in 1:p) {
      betastar[i] <- t(X[,i])%*%(y-X[,-i]%*%beta[-i])/crossprod(X[,i])
      if(1/crossprod(X[,i]) *N*lambda/2 > abs(betastar[i])){
        beta[i] = 0
      }else if(betastar[i] > 0){
        beta[i] = betastar[i] - 1/crossprod(X[,i])*N*lambda/2
      }else{
        beta[i] = betastar[i] + 1/crossprod(X[,i])*N*lambda/2
      }
      
    }
    ###Stopping rule
    if(sum(((beta0old-beta)/beta)^2, na.rm = T)<1e-5){
      break;
    }
  }

###Now we fit LASSO using glmnet.
###As Yulin mentioned the solution will match with glmnet for lambda/2 due to 
###the nature of the glmnet LOSS function (they have an extra 1/2 in the first term https://glmnet.stanford.edu/articles/glmnet.html)

lambda = lambda/2
fit <- glmnet::glmnet(X, y, alpha = 1, intercept = F, lambda = lambda)
betalasso <- as.matrix(fit$beta)

###Checking closeness between our own implementation and glmnet
range(beta-betalasso)
