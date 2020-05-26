# p <- 1
# q <- 1
# n <- 10000
# fi <- c(0.5)
# teta <- c(-0.17)
# 
# 
# #autoregression of moving average process model
# y <- rep(0, n)
# eps <- rnorm(n, mean = 0, sd = 5^(1/2))
# for (t in (p+1):n){
#   y[t] <- sum(fi*y[(t-1):(t-p)]) + eps[t] + sum(-teta*eps[(t-1):(t-q)])
# }
# 
# # d <- 2
# x <- numeric(length = n-2)
# for (i in 3:n)
#   x[i] <- y[i]+2*x[i-1]-x[i-2]
# n <- n-2
# =============================================
x <- read.csv("z3.csv")$x
n <- length(x)
p <- 1
q <- 1


# # d <- 1
# y <- x[2:n]-x[1:(n-1)]
# n <- n-1
# 
# autocov <- vector(length = n)
# for (k in 1:n)
#   autocov[k] <- mean(y[1:(n-k+1)]*y[k:n])
# 
# autocor <- autocov/autocov[1]
# 
# plot(abs(autocor[1:10]), type = "l", ylim = c(0,1))

#d <- 2
y <- x[3:n]-2*x[2:(n-1)]+x[1:(n-2)]
n <- n-2

autocov <- vector(length = n)
for (k in 1:n)
  autocov[k] <- mean(y[1:(n-k+1)]*y[k:n])

autocor <- autocov/autocov[1]

plot(abs(autocor[1:10]), type = "l", ylim = c(0,1))

# #d <- 3
# y <- x[4:n]-3*x[3:(n-1)]-x[2:(n-2)]-x[1:(n-3)]
# n <- n-3
# 
# autocov <- vector(length = n)
# for (k in 1:n)
#   autocov[k] <- mean(y[1:(n-k+1)]*y[k:n])
# 
# autocor <- autocov/autocov[1]
# 
# plot(abs(autocor[1:10]), type = "l", ylim = c(0,1))

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
  w[t] <- y[t]-sum(fiest*y[(t-1):(t-p)])

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

eps_est <- rep(0, n+l)
for (t in (q+1):n)
  eps_est[t] <- y[t]+sum(tetaest*eps_est[(t-1):(t-q)])-sum(fiest*y[(t-1):(t-p)])
l <- 5
eps_est[(n+1):(n+l)] <- 0

y <- c(y, rep(0,l))
for (i in 1:l){
  y[n+i] <- sum(fiest*y[(n+i-1):(n+i-p)])-sum(tetaest*eps_est[(n+i-1):(n+i-q)])+eps_est[n+i]
}


# d <- 2
n <- n+2
x_est <- vector(length = n+l)
x_est[1:n] <- x[1:n]
for (i in 1:l)
  x_est[n+i] <- y[n+i-2]+2*x_est[n+i-1]-x_est[n+i-2]


# d <- 1
n <- n+1
x_est <- vector(length = n+l)
x_est[1:n] <- x[1:n]
# x_est <- x[1:n]
# for (i in 1:l)
#   x_est[n+i] <- y[n+i-1]+x_est[n+i-1]


# d <- 3
n <- n+3
x_est <- vector(length = n+l)
x_est[1:n] <- x[1:n]
# x_est <- x[1:n]
# for (i in 1:l)
#   x_est[n+i] <- y[n+i-3]+3*x_est[n+i-1]+x_est[n+i-2]+x_est[n+i-3]
