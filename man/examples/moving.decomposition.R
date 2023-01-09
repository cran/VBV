### Usage of moving.decomposition

t <- 1:121 # equidistant time points, i.e. 5 days

m <- 11

p <- 2     # maximally quadratic
q <- c(1, 3, 5)   # 'seasonal' components within the base period
base.period <- 24 # i.e. hourly data with daily cycles
l1 <- 1    
l2 <- 1

m.dec <- moving.decomposition( length(t), p, q, m, base.period, l1, l2)
