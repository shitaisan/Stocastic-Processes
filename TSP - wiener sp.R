tn <- 0.5
h <- 0.00001
n <- 50
x <- replicate(n, cumsum(rnorm(tn/h, sd = h^(1/2)))) # генерируем n траекторий ВП 
print(ks.test(x[tn/h,], pnorm, sd = tn^(1/2)))

a <- 0.01
# функция для определния момента первого пересечения процессов уровня а
first <- function(x){
  for (i in 1:length(x))
    if (x[i]>=a)
      return (i)
  # если не вышли внутри for - вся траектория х<a, продолжаем ее, пока не x>=a
  while (x[length(x)]<a)
    x <- c(x, x[length(x)]+rnorm(1, sd = h^(1/2)))
  return (length(x))
}

# ф.р. для тау
distrtau <- function(t){
  return (2-2*pnorm(a/t^(1/2)))
}

# каждый столбец матрицы х - траектория ВП, применяем first() для каждой траектории
# моменты времени = index*h
tau <- apply(x, 2, first)*h
print(ks.test(tau, distrtau))

