### Usage of decomposition
t <- 1:121 # equidistant time points, i.e. 5 days
p <- 2     # maximally quadratic
q <- c(1, 3, 5)   # 'seasonal' components within the base period
base.period <- 24 # i.e. hourly data with daily cycles
l1 <- 1    
l2 <- 10

dec <- decomposition( t, p, q, base.period, l1, l2)
### Note: decomosition is independent of data, only depends on time
