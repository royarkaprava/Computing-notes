PreLinearCG <- function(minv, A, b, x0, tol=1e-5){
  xk = x0;
  rk = b - A %*% xk;
  zk = minv * rk;
  dk = zk;
  error = sum(rk*rk);
  
  adk = A %*% dk;
  rkzk = sum(rk*zk);
  
  alpha = rkzk / sum(dk*adk);
  
  beta = sum(rk*rk) / rkzk;
  
  n = nrow(A);
  
  num_iter = 0;
  alphadkvec <- 0
  while(error > tol){
    adk = A %*% dk;
    rkzk = sum(rk*zk);
    
    alpha = rkzk / sum(dk*adk);
    xk = xk + alpha * dk;
    rk = rk - alpha * adk;
    zk = minv * rk;
    
    beta = sum(rk*zk) / rkzk;
    dk = zk + beta * dk;
    
    num_iter = num_iter + 1;
    alphadkvec[num_iter] <- alpha*sqrt(sum(dk^2))
    if(num_iter > n){break;}
    error = sum(rk*rk);
  }
  
  return(list(sol=xk, alpha=alphadkvec))
}
