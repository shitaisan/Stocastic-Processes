n <-100

# autoregression process model
p <- 2
fi <- c(0.5, 0.2)
AR <- rep(0, n)
eps <- rnorm(n, mean = 0, sd = (2.36)^(1/2))
for (t in (p+1):n){
  AR[t] <- eps[t]+sum(fi[1:p]*AR[(t-1):(t-p)])
}

# moving average process model
q <- 2
teta <- c(-0.4, 0.05)
MA <- rep(0, n)
eps <- rnorm(n, mean = 0, sd = 1)
for (i in (q+1):n){
  MA[i] <- eps[i]-sum(teta*eps[(i-1):(i-q)])
}

#autoregression of moving average process model
ARMA <- rep(0, n)
eps <- rnorm(n, mean = 0, sd = 1)
for (t in (p+1):n){
  ARMA[t] <- sum(fi*ARMA[(t-1):(t-p)]) + eps[t] + sum(-teta*eps[(t-1):(t-q)])
}

x <- ARMA
autocov <- vector(length = n)
for (k in 1:n)
  autocov[k] <- mean(x[1:(n-k+1)]*x[k:n])

autocor <- autocov/autocov[1]

plot(abs(autocor[1:10]), type = "l")

privautocor <- vector(length = n)
# for (k in 1:(n-1)){
#   A <- vector(length=k*k)
#   dim(A) <- c(k, k)
#   for (i in 1:k)
#     for (j in 1:k)
#       A[i,j] <- autocor[abs(i-j)+1]
#   B <- A
#   B[, k] <- autocor[2:(k+1)]
#   privautocor[k] <- det(B)/det(A)
# }
for (k in 1:10){
  A <- vector(length=k*k)
  dim(A) <- c(k,k)
  for (i in 1:k)
    for (j in 1:k)
      A[i,j] <- autocor[abs(i-j)+1]
  b <- vector(length=k)
  b <- autocor[2:(k+1)]
  privautocor[k] <- solve(a = A, b = b)[k] 
}

lines(abs(privautocor[1:10]), col = 'blue')
sdCor <- vector(length = n)
for (k in 1:n)
  sdCor[k] <- (1/n*(1+2*(sum(autocor[2:k]^2))))^(1/2)
sdPrCor <- rep((1/n)^(1/2), k)

lines(sdCor, lty = 3, col = 'red')
lines(sdPrCor, lty = 3, col = 'green')
  
  
  
  
  
  
  
  















