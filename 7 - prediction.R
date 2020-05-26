x <- read.csv("3.csv", sep = ',', dec = '.')$x
n <- length(x)
t <- 1:n
tr <- lm(x~t)

# h <- seq(0, 1, by = 0.001)
# for (i in 1:length(h)){
#   a[i] <- 2/n*sum(x*cos(2*pi*t*h[i]))
#   b[i] <- 2/n*sum(x*sin(2*pi*t*h[i]))
# }
# p <- a^2+b^2
# plot(p)
# J <- 5
# 
# seas <- rep(0, n)
# for (j in J)
#   seas <- seas+a[j]*cos(2*pi*t*h[j])+b[j]*sin(2*pi*t*h[j])


h <- 1/300
a <- 2/n*sum(x*cos(2*pi*t*h))
b <- 2/n*sum(x*sin(2*pi*t*h))
seas <- a*cos(2*pi*t*h)+b*sin(2*pi*t*h)
y <- x - tr$fitted.values - seas
plot(y)

autocov <- vector(length = n)
for (k in 1:n)
  autocov[k] <- mean(y[1:(n-k+1)]*y[k:n])

autocor <- autocov/autocov[1]
plot(abs(autocor[1:10]), type = "l", ylim = c(0,1))

privautocor <- vector(length = 10)
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

lines(abs(privautocor), col = 'blue')

sdCor <- vector(length = n)
for (k in 1:n)
  sdCor[k] <- (1/n*(1+2*(sum(autocor[2:k]^2))))^(1/2)
sdPrCor <- rep((1/n)^(1/2), k)
lines(sdCor, lty = 3, col = 'red')
lines(sdPrCor, lty = 3, col = 'green')


legend(6, 1, legend = c("autocor", "priv autocor", "sd for autocor", "sd for priv autocor"), col = c(1, "blue", "red", "green"), lty = c(1, 1, 3, 3))


p <- 1
q <- 1

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
acc <- 10^(-3)
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

l <- 5
eps_est <- rep(0, n+l)
for (t in (q+1):n)
  eps_est[t] <- y[t]-sum(fiest*y[(t-1):(t-p)])+sum(tetaest*eps_est[(t-1):(t-q)])


y <- c(y, rep(0,l))
for (i in 1:l){
  y[n+i] <- sum(fiest*y[(n+i-1):(n+i-p)])-sum(tetaest*eps_est[(n+i-1):(n+i-q)])+eps_est[n+i]
}


t <- 1001:1005
newseas <- a*cos(2*pi*t*h)+b*sin(2*pi*t*h)
newtr <- tr$coefficients[1]+t*tr$coefficients[2]
x <- c(x, newseas+newtr+y[t])


  