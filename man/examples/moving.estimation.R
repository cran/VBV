### usage of moving.estimation

t <- 1:121 # equidistant time points, i.e. 5 days
y <- 0.1*t + 0.02*t^2 + rnorm(length(t))

m <- 11

p <- 2     # maximally quadratic
q <- c(1, 3, 5)   # 'seasonal' components within the base period
base.period <- 24 # i.e. hourly data with daily cycles
l1 <- 1    
l2 <- 1

m.est <- moving.estimation( t, y, p, q, m, base.period, l1, l2)
plot(m.est$data)
lines(m.est$trend + m.est$season)

