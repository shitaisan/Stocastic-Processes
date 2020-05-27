q <- 3
n <- 10000
teta <- c(-0.7, -0.5, 0.01) # sum of abs < 1


acc <- 10^(-9)

#moving average process model
x <- rep(0, n)
eps <- rnorm(n, mean = 0, sd = 2)
for (i in (q+1):n){
  x[i] <- eps[i]-sum(teta*eps[(i-1):(i-q)])
}
#x <- as.vector(read.csv('ma.csv')$x)


#Moving average parameters teta estimation
autocov <- vector(length = q+1)
for (k in 0:q)
  autocov[k+1] <- mean(x[1:(n-k)]*x[(k+1):n])

tetaest <- rep(0, q)
varest <- 1
varestPrev <- acc
tetaestPrev <- rep(acc, q)
while (abs(varestPrev-varest)>=acc || max(abs(tetaest-tetaestPrev))>=acc) {
  varestPrev <- varest
  tetaestPrev <- tetaest
  varest <- autocov[1]/(1+sum(tetaest^2))
  tetaest[q] <- -autocov[q+1]/varest
  tetaest[q-1] <- -autocov[q]/varest + tetaest[1]*tetaest[q]
  if (q>=3){
    for (i in (q-2):1)
      tetaest[i] <- -autocov[i+1]/varest + sum(tetaest[1:(q-i)]*tetaest[(i+1):q])
  }
}
