p <- 2
n <- 10000

fi <- c(0.348, -0.626)


#autoregression process model
x <- rep(0, n)
eps <- rnorm(n, mean = 0, sd = (2.36)^(1/2))
for (t in (p+1):n){
  x[t] <- eps[t]+sum(fi[1:p]*x[(t-1):(t-p)])
}
#x <- as.vector(read.csv('ar.csv')$x)


#autocovariance estimations
autocov <- vector(length = p+1)
for (k in 0:p)
  autocov[k+1] <- mean(x[1:(n-k)]*x[(k+1):n])


#autoregression parameters estimations from Ax=b solution
A <- vector(length=p*p)
dim(A) <- c(p,p)
for (i in 1:p)
  for (j in 1:p)
    A[i,j] <- autocov[abs(i-j)+1]
b <- vector(length=p)
b <- autocov[2:(p+1)]

fiest <- solve(a = A, b = b) 
varest <- autocov[1]-sum(fiest*autocov[2:(p+1)])

