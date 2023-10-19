LinearCG <- function(A, b, x0, tol=1e-5){
  xk = x0
  rk = b - A %*% xk
  dk = rk
  error = sum(rk^2)
  
  n <- nrow(A)
  
  num_iter = 0
  alphadkvec <- 0
  while(error > tol){
    adk = A %*% dk
    rk2old = sum(rk^2)
    
    alpha = rk2old / sum(dk*adk)
    
    xk = xk + alpha * dk
    rk = rk - alpha * adk
    beta = sum(rk^2) / rk2old
    dk = rk + beta * dk
    
    num_iter = num_iter + 1
    alphadkvec[num_iter] <- alpha*sqrt(sum(dk^2))
    if(num_iter > n){break}
    error <- sum(rk^2)
  }
  
  rk_norm = sqrt(sum(rk^2))
  
  return(list(sol=xk, alpha=alphadkvec))
}
Rcpp::sourceCpp("Conjugate_gradient.cpp")

p <- 1000
n <- 1500
X = matrix(rnorm(n*p), n, p)
beta <- rnorm(p)
y = X%*%beta + rnorm(n, 0, 1)


A <- crossprod(X);b <- crossprod(X, y); 
t1 <- proc.time()
system.time(x3 <- LinearCGc(A, b, rep(0, p), tol=1e-4))
t2 <- proc.time()
mean((x3-beta)^2)

t1 <- proc.time()
x2 <- LinearCG(A, b, rep(0, p), tol=1e-4)
t2 <- proc.time()
t2-t1
plot(x2$alpha)
mean((x2$sol-beta)^2)

system.time(x1 <- lm(y~X-1)$coefficients)
mean((x1-beta)^2)

system.time(x0 <- solve(A, b))
mean((x0-beta)^2)
