B <- 1000
nvec <- c(10, 20, 30, 50, 100,500, 1000,2000)
dat_list <- list()
for (i in 1:length(nvec)) {
  n <- nvec[i]
  # alpha = 1, beta = 2
  data <- matrix(rgamma(B*n, shape = 1, rate = 2), nrow = B, ncol = n)
  Sn <- rowMeans(data)
  
  plot(density(Sn), main = paste("Density of Sn, n = ", n))
  #points(density(normal_sample), type='l', col=2)
  
  qqnorm(Sn, pch = 1, frame = FALSE)
  qqline(Sn, col = "steelblue", lwd = 2)
  
  print(shapiro.test(Sn)$p.value)
}


B <- 1000
nvec <- c(10, 20, 30, 50, 100, 150, 200, 500, 1000,2000)

power_table <- data.frame(matrix(ncol = 3, nrow = length(nvec)))
for(i in 1:length(nvec)){
  n <- nvec[i]
  xnull <- rgamma(n*B, shape = 270, rate = 2)
  xnull <- matrix(xnull, nrow = B, ncol = n)
  Sn_null <- rowMeans(xnull)
  type1_error <- mean(abs(Sn_null - 135) > 2)
  xalt <- rgamma(n*B, shape = 275, rate = 2)
  xalt <- matrix(xalt, nrow = B, ncol = n)
  Sn_alt <- rowMeans(xalt)
  power <- mean(abs(Sn_alt - 135) > 2)
  power_table[i,] <- c(n, type1_error, power)
}
colnames(power_table) <- c('n', 'Type-I error', 'power')

power_table
