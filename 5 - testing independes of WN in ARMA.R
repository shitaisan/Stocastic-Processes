acc <- 10^(-9)


x <- read.csv("3.csv")$x #if x is MA(2)
n <- length(x)
q <- 2

autocov <- vector(length = q+1)
for (k in 0:q)
  autocov[k+1] <- mean(x[1:(n-k)]*x[(k+1):n])

autocor <- autocov/autocov[1]

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

eps_est <- rep(0, n)
for (t in (q+1):n)
  eps_est[t] <- x[t]+sum(tetaest*eps_est[(t-1):(t-q)])

autocov <- vector(length = n)
for (k in 1:n)
  autocov[k] <- mean(eps_est[1:(n-k+1)]*eps_est[k:n])

autocor <- autocov/autocov[1]

K <- 4
# testing wn independence with Box-Piarson statistic
alpha <- 0.05
Q <- n*sum(autocor[2:K]^2)
pval <- 1-pchisq(Q, K-q)
if (pval>alpha){
  cat("can't reject independence, pval = ", pval, "is too much")
} else {cat("reject null with pval = ", pval)}
