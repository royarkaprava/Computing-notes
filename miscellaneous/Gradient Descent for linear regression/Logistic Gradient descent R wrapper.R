Rcpp::sourceCpp("Logistic_regression_gradient_descent.cpp")

x = rnorm(1000, 0, 1)
y = rbinom(1000, 1, 1/(1+exp(3+2*x)))

X <- cbind(rep(1, 1000), x)

out = logisticGDfit(X, y, 10000, 0.1, 0)

out = logisticGDfitR(X, y, 10000, 0.1, 0)

out$beta
plot(out$loss[which(out$loss>0)], type='l')

mylogit <- glm(y ~ X-1, family = "binomial")
mylogit$coefficients #The coefficient are negative to what we got as glm() sets pi=1/(1+exp(-Xbeta))

#As a practice change L (the fourth entry) form 0.1 and Total_itr form 10000 to see their effects.


#Logistic regression

mydata <- as.matrix(read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv", header=T))
y <- mydata[, 1]
X <- mydata[, -1]

mylogit <- glm(y ~ scale(X)-1, family = "binomial")
mylogit$coefficients #The coefficient are negative to what we got as glm() sets pi=1/(1+exp(-Xbeta))

out = logisticGDfit(scale(X), y, 1000, 0.1, 0)

out$beta