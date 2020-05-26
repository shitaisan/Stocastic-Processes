df <-  read.csv("z1-1.csv")
plot(df$x, df$y, type="l")
x <- df$x
y <- df$y
rng <- function(ii, maxi){
  return (ii[ii>=1 & ii<=maxi])
}

delta <- 2
z1 <- c()
for (i in 1:length(y)){
  z1[i] <- mean((y[rng((i-delta):(i+delta), length(y))]))
}
lines(x, z1, col = "blue")

z2 <- c()
for (i in 1:length(y)){
  z2[i] <- median((y[rng((i-delta):(i+delta), length(y))]))
}
lines(x, z2, col = "green")

df2 <- read.csv("z1-2.csv")
tr <- lm(df2$y~df2$x)$fitted.values
plot(df2, type='l')
lines(df2$x, tr, col = 'blue')
eps <- df2$y-tr
lines(df2$x, eps, type='l', col = "red")
mean(eps^2)
