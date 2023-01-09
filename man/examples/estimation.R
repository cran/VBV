### using of estimation

t <- 1:121 # equidistant time points, i.e. 5 days
y <- 0.1*t + sin(t) + rnorm(length(t))

p <- 2     # maximally quadratic
q <- c(1, 3, 5)   # 'seasonal' components within the base period
base.period <- 24 # i.e. hourly data with daily cycles
l1 <- 1    
l2 <- 10

est <- estimation( t, y, p, q, base.period, l1, l2)
plot(est$data)
lines(est$trend + est$season)
