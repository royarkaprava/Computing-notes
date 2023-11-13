#Slightly modified the 

n = 10
a = 2
b = 3
# Full Gibbs
X <- 0 ; Y <- 0
X[1] <- 1 #init value for x_0
Total_itr <- 2000
set.seed(123456)
for(i in 1:Total_itr){
  # sampling from Y|X
  Y[i] <- rbeta(1,X[i]+a,n-X[i]+b)
  # sampling from X|Y
  X[i+1] <- rbinom(1,n,Y[i])
}
out <- cbind(x=X[-1],y=Y)
outp <- out[-c(1:500), ]

colMeans(outp)
##########################MH sampling with transformation#######################

set.seed(123456)
fderiv = function(l){
  deriv = exp(l)/(1+exp(l))^2
  return (deriv)
}
finv = function(lambda){
  y = exp(lambda)/(1+exp(lambda))
  return (y)
}
ar1 <- 0; sd1 = 1
XMH <- 1; lambdas = 0; YMH=finv(0)
a <- 2; b <- 3

Total_itr <- 10000
for(i in 1:Total_itr){
  xold = XMH[i]
  lambdaold = lambdas[i]
  lambdac <- lambdaold + rnorm(n=1,sd=sd1)
  yold = finv(lambdaold)
  yc = finv(lambdac)
  # likelihood for lambda with Jacobina adjustment
  lold = (xold+a-1)*log(yold)+(n-xold+b-1)*log(1-yold)+log(fderiv(lambdaold))
  lnew = (xold+a-1)*log(yc)+(n-xold+b-1)*log(1-yc)+log(fderiv(lambdac))
  R <- lnew - lold
  if(is.na(R)) {R = -Inf}
  if(is.nan(R)) {R = -Inf}
  logu <- log(runif(1))
  if (logu < R){
    lambdas[i+1] <- lambdac # new sample of lambda
    YMH[i+1] = yc #new sample of Y is yc
    ar1 <- ar1 + 1
  }
  if (logu >= R) {
    lambdas[i+1] <- lambdaold
    YMH[i+1] = yold
  }
  XMH[i+1] <- rbinom(1,n,YMH[i+1]) #update x
  if(i %% 100 == 0){
    #Adjust the standard deviation to maintain an acceptance rate within (0.2,0.5)
    ar <- ar1 / i
    if(ar>.50){sd1 <- sd1 * 2}
    if(ar<.20){sd1 <- sd1 / 2}
  }
}
outMH <- cbind(x=XMH,y=YMH)
outpMH <- outMH[-c(1:5000), ]
acf(outpMH[,1],lag.max = 100)

colMeans(outpMH)


##########################MH sampling without transformation#######################

set.seed(123456)

ar1 <- 0; sd1 = 0.1
XMH <- 1; lambdas = 0; YMH=0.5
a <- 2; b <- 3

Total_itr <- 10000
for(i in 1:Total_itr){
  xold = XMH[i]
  yold = YMH[i]
  yc = yold + rnorm(n=1,sd=sd1)
  # likelihood for lambda with Jacobina adjustment
  lold = (xold+a-1)*log(yold)+(n-xold+b-1)*log(1-yold)+log(fderiv(lambdaold))
  lnew = (xold+a-1)*log(yc)+(n-xold+b-1)*log(1-yc)+log(fderiv(lambdac))
  R <- lnew - lold
  
  #If yc falls outside of (0,1), lnew will automatically be NA and the sample will not be accepted
  
  if(is.na(R)) {R = -Inf}
  if(is.nan(R)) {R = -Inf}
  logu <- log(runif(1))
  if (logu < R){
    YMH[i+1] = yc #new sample of Y is yc
    ar1 <- ar1 + 1
  }
  if (logu >= R) {
    YMH[i+1] = yold
  }
  XMH[i+1] <- rbinom(1,n,YMH[i+1]) #update x
  if(i %% 100 == 0){
    #Adjust the standard deviation to maintain an acceptance rate within (0.2,0.5)
    ar <- ar1 / i
    if(ar>.50){sd1 <- sd1 * 2}
    if(ar<.20){sd1 <- sd1 / 2}
  }
}
outMH <- cbind(x=XMH,y=YMH)
outpMH <- outMH[-c(1:5000), ]
acf(outpMH[,1],lag.max = 100)

colMeans(outpMH)