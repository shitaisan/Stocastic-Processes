p <- 1
q <- 1
n <- 10000
fi <- c(0.5)
teta <- c(-0.17)


#autoregression of moving average process model
x <- rep(0, n)
eps <- rnorm(n, mean = 0, sd = 5^(1/2))
for (t in (p+1):n){
  x[t] <- sum(fi*x[(t-1):(t-p)]) + eps[t] + sum(-teta*eps[(t-1):(t-q)])
}
#x <- as.vector(read.csv('ma.csv')$x)

#autocovariance estimations
autocov <- vector(length = p + q + 1)
for (k in 0:(p + q))
  autocov[k+1] <- mean(x[1:(n-k)]*x[(k+1):n])

#autoregression parameters estimations from Ax=b solution
A <- vector(length = p*p)
dim(A) <- c(p,p)
for (i in 1:p)
  for (j in 1:p)
    A[i,j] <- autocov[q+abs(i-j)+1]
b <- vector(length=p)
b <- autocov[(q+2):(q+p+1)]

fiest <- solve(a = A, b = b) 

w <- rep(0, n)
#removing autoregression from our ARMS
for (t in (p+1):n)
  w[t] <- x[t]-sum(fiest*x[(t-1):(t-p)])

#estimation for autocovariance w(t)
wautocov <- vector(length = q+1)
for (k in 0:q)
  wautocov[k+1] <- mean(w[1:(n-k)]*w[(k+1):n])

#estimation for w(t) parameters teta and variance as for moving average
acc <- 10^(-9)
tetaest <- rep(0, q)
varest <- 1
varestPrev <- acc
tetaestPrev <- rep(acc, q)
while (abs(varestPrev-varest)>=acc || max(abs(tetaest-tetaestPrev))>=acc) {
  varestPrev <- varest
  tetaestPrev <- tetaest
  varest <- wautocov[1]/(1+sum(tetaest^2))
  tetaest[q] <- -wautocov[q+1]/varest
  if (q > 1) 
    tetaest[q-1] <- -wautocov[q]/varest + tetaest[1]*tetaest[q]
  if (q > 2){
    for (i in (q-2):1)
      tetaest[i] <- -wautocov[i+1]/varest + sum(tetaest[1:(q-i)]*tetaest[(i+1):q])
  }
}
                   
