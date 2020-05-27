df <- read.csv("daily-min-temperatures.csv", sep = ",", dec = ".")
a <- c()
b <- c()
n <- nrow(df)
for (j in 1:n){
  a[j] <- 2/n*sum(df$Temp*cos(2*pi*df$Date*j/n))
  b[j] <- 2/n*sum(df$Temp*sin(2*pi*df$Date*j/n))
}
p <- a^2+b^2
plot(p[(n-10):n])
plot(df)
j <- 10
seas <- c()
for (t in 1:n){
  seas[t] <- a[j]*cos(2*pi*j/n*df[1,t])+b[j]*sin(2*pi*j/n*df[1,t])
}
