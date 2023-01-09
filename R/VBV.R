#' decomposition - decompose a time series with VBV
#'
#' @param t.vec vector of observation points.
#' @param p maximum exponent in polynomial for trend
#' @param q.vec vector containing frequencies to use for seasonal component,
#' given as integers, i.e. c(1, 3, 5) for 1/2*pi, 3/2*pi, 5/2*pi (times length of base period)
#' @param base.period base period in number of observations, i.e. 12 for monthly data with yearly oscillations
#' @param lambda1 penalty weight for smoothness of trend
#' @param lambda2 penalty weight for smoothness of seasonal component
#' (lambda1 == lambda2 == Inf result in estimations of the original Berliner Verfahren)
#' @return list with the following components:
#' \itemize{
#' \item{trend}{A function which returns the appropriate weights if applied to a point in time}
#' \item{saison}{A function which returns the appropriate weights if applied to a point in time}
#' \item{A, G1, G2}{Some matrices that allow to calclate SSE etc. Exposed only to reuse their calculation. See the referenced paper for details.}
#' }
#' @example man/examples/decomposition.R
#' @export

decomposition <- function(t.vec, p , q.vec, base.period, lambda1, lambda2) {
    mitte <- sum(range(t.vec))/2
    t.vec <- t.vec - mitte
    n <- length(t.vec)

    f10  <- function(t, p, q){
        c(t^(0:(p-1)), rep(0,2*q) )
    }

    f02  <- function(t, p, q.vec, base.period){
        q <- length(q.vec)
        erg  <- rep(0,p+2*q)
        erg[ seq( p+1, p+2*q-1, 2)]  <- cos(2*pi*t*q.vec/base.period)
        erg[ seq( p+2, p+2*q, 2)]  <- sin(2*pi*t*q.vec/base.period)
        erg    
    }

    g1.plus  <- function(tk, t, p){
        if (t <= tk) return (0)
        (t-tk)^(2*p -1)
    }
    
    g1 <- function(t, t.vec, p){
        sapply(t.vec, g1.plus, t=t, p=p)
    }
    
    g2.plus <- function(tk, t, base.period, q.vec){
        if (t <= tk) return(0)

        r.vec <- q.vec/mean(q.vec)
        r2.vec  <- r.vec^2
        
        d.vec  <- 1/r.vec
        c.vec  <- 1/r.vec
        
        q <- length(q.vec)
        
        for (i in 1:q)
        {
            d.vec[i]  <- d.vec[i] - 4*r.vec[i] * sum(1/(r2.vec[-i]-r2.vec[i]))
            c.vec[i]  <- c.vec[i] / prod(r2.vec[-i]-r2.vec[i])
        }
        
        sum(c.vec^2*( d.vec * sin( 2*pi* (t-tk)*q.vec/base.period) -
                              mean(q.vec) * 2 * pi * (t -tk)/base.period *
                              cos(2*pi*(t-tk)*q.vec/base.period)  )
            )
    }

    g2  <- function(t, t.vec, q.vec, base.period){
        d <- 1/(2*(2*pi*mean(q.vec)/base.period)^(4*length(q.vec)-1))
        d*sapply(t.vec, g2.plus, t=t, q.vec=q.vec, base.period=base.period)
    }

    FGES  <- t(
               sapply(t.vec, f10, p=p, q=length(q.vec)) +
               sapply(t.vec, f02, p=p, q.vec=q.vec, base.period=base.period)
               )
     
    Bstar <- qr.solve( (t(FGES) %*% FGES), tol=1e-17) %*% t(FGES)
    Astar <- diag(1, length(t.vec)) - FGES %*% Bstar

    G1 <-  sapply(t.vec, g1, t.vec=t.vec, p=p)
    G2 <-  sapply(t.vec, g2, t.vec=t.vec, q.vec=q.vec, base.period=base.period)
    GGES  <- -t( G1/lambda1 + G2/lambda2 )

    A <- qr.solve(diag(1, n) - Astar %*% GGES, tol=1e-17) %*% Astar
    B <- Bstar%*% (diag(1, n)+ GGES %*% A)
    
    list(
        trend = function(t) {
            f10(t,p,length(q.vec))%*% B  +
                g1(t, t.vec,p)%*%((1/lambda1)*A)},
        season = function(t) {
        f02(t,p,q.vec,base.period) %*%  B +
                g2(t, t.vec, q.vec, base.period) %*% ((1/lambda2)*A)},
        A = A , Astar = Astar, G1 = G1, G2 = G2 )
}

#' estimation -- estimate trend and seasonal components statically
#'@param t.vec vector of points in time as integers
#'@param y.vec vector of data
#'@param p maximum exponent in polynomial for trend
#'@param q.vec vector containing frequencies to use for seasonal component, given as integers, i.e. c(1, 3, 5) for 1/2*pi, 3/2*pi, 5/2*pi (times length of base period)
#'@param base.period base period in number of observations, i.e. 12 for monthly data with yearly oscillations
#'@param lambda1 penalty weight for smoothness of trend
#'@param lambda2 penalty weight for smoothness of seasonal component (lambda1 == lambda2 == Inf result in estimations of the original Berliner Verfahren)
#'@return A dataframe with the following components:
#' \itemize{
#' \item{data}{original data y.vec}
#' \item{trend}{vector of estimated trend of length length(y.vec)}
#' \item{season}{vector of estimated season of length length(y.vec)}
#' }
#' @example man/examples/estimation.R
#' @export
estimation <- function(t.vec, y.vec, p, q.vec,  base.period, lambda1, lambda2){
    t.vec <- t.vec - sum(range(t.vec))/2

    n <- length(t.vec)
    trend  <- rep(0, n)
    season <- rep(0, n)

    parts <- decomposition(t.vec, p, q.vec, base.period, lambda1, lambda2)

    for (t in 1:n) {
        trend[t]  <- parts$trend(t.vec[t]) %*% y.vec
        season[t] <- parts$season(t.vec[t]) %*% y.vec
    }

    data.frame(data=y.vec, trend = trend, season = season)
}

#'moving.estimation -- estimate locally optimized trend and season figures
#'@param t.vec vector of points in time as integers
#'@param y.vec vector of data
#'@param p maximum exponent in polynomial for trend
#'@param q.vec vector containing frequencies to use for seasonal component, given as integers, i.e. c(1, 3, 5) for 1/2*pi, 3/2*pi, 5/2*pi (times length of base period)
#'@param m width of moving window
#'@param base.period base period in number of observations, i.e. 12 for monthly data with yearly oscillations
#'@param lambda1 penalty weight for smoothness of trend
#'@param lambda2 penalty weight for smoothness of seasonal component
#'@note lambda1 == lambda2 == Inf result in estimations of the original Berliner Verfahren
#'@return A dataframe with the following components:
#'\itemize{
#' \item{data}{original data y.vec}
#' \item{trend}{vector of estimated trend of length length(y.vec)}
#' \item{season}{vector of estimated season of length length(y.vec)}
#' }
#' @export
moving.estimation <- function( t.vec, y.vec, p, q.vec, m, base.period, lambda1, lambda2){
    weights <- moving.decomposition( length(t.vec), p, q.vec, m, base.period, lambda1, lambda2)
    trend <-  weights$W1 %*% y.vec
    dimnames(trend) <- NULL
    season <- weights$W2 %*% y.vec
    dimnames(season) <- NULL
    dimnames(y.vec) <- NULL
    data.frame(data=y.vec, trend = trend, season = season)
}

#' moving.decomposition -- decompose a times series into locally estimated trend and season figures
#'
#' @param n number of observation points (must be odd!). Internally this will be transformed to seq( -(n-1)/2, (n-1)/2, 1)
#' @param p maximum exponent in polynomial for trend
#' @param q.vec vector containing frequencies to use for seasonal component, given as integers, i.e. c(1, 3, 5) for 1/2*pi, 3/2*pi, 5/2*pi (times length of base period)
#' @param m width of moving window
#'@param base.period base period in number of observations, i.e. 12 for monthly data with yearly oscillations
#'@param lambda1 penalty weight for smoothness of trend
#'@param lambda2 penalty weight for smoothness of seasonal component
#'@note lambda1 == lambda2 == Inf result in estimations of the original Berliner Verfahren
#' @return list with the following components:
#'\itemize{
#'   \item{W1}{nxn matrix of weights. Trend is estimated as W1 %*% y, if y is the data vector}
#'   \item{W2}{nxn matrix of weights. Season is estimated as W2 %*% y, if y is the data vector}
#' }
#' @example man/examples/moving.decomposition.R
#' @export
moving.decomposition <- function(n, p, q.vec, m, base.period, lambda1, lambda2){
    if (!(n%%2==1)) {
       warning("n must be odd!")
       stop()
    }
    k <- (m-1)/2
    W1 <- matrix(0,nrow=n, ncol=n)
    W2 <- matrix(0,nrow=n, ncol=n)

    parts <- decomposition(1:m, p, q.vec, base.period, lambda1, lambda2)

    for ( zeile in -k:0) {
        W1[zeile + 1 + k, 1:m] <- parts$trend(zeile)
        W2[zeile + 1 + k, 1:m] <- parts$season(zeile)
    }

    zeile1 <- W1[k+1,1:m]
    zeile2 <- W2[k+1,1:m] 

    for (zeile in (k + 2):(n - k -1 )) {
        W1[zeile, (zeile - k ):( zeile- k -1 +m)  ] <- zeile1
        W2[zeile, (zeile - k ):( zeile- k -1 +m)  ] <- zeile2
    }

    for ( zeile in (n-k):n ) {
        W1[zeile, (n-m+1):n] <- parts$trend(zeile-n+k)
        W2[zeile, (n-m+1):n] <- parts$season(zeile-n+k)
    }
    
    list ( W1=W1, W2=W2)
}

