
# Expectations ----------------------------------------------------------------------------------------------------

## see Appendix A1
EZ0ab <- function(a, b, t, beta){
    ##
    stopifnot(beta<0)
    stopifnot(t>0)
    stopifnot(a<=b || b==Inf)
    stopifnot(a>=0)
    if(a==b) return(0)
    ##
    p1 <- if (b==Inf) Inf else (0.5*beta^2*t+log(b))/(-beta*sqrt(t))
    p2 <- if (a==0) -Inf else if (a==Inf) Inf else (0.5*beta^2*t+log(a))/(-beta*sqrt(t))
    ##
    pnorm(p1) - pnorm(p2)
}
EZ1ab <- function(a, b, t, beta){
    ##
    stopifnot(beta<0)
    stopifnot(t>0)
    stopifnot(a<=b || b==Inf)
    stopifnot(a>=0)
    if(a==b) return(0)
    ##
    p1 <- if (b==Inf) Inf else (0.5*beta^2*t-log(b))/(beta*sqrt(t))
    p2 <- if (a==0) -Inf else if (a==Inf) Inf else (0.5*beta^2*t-log(a))/(beta*sqrt(t))
    ##
    pnorm(p1) - pnorm(p2)
}
EZ2ab <- function(a, b, t, beta){
    ##
    stopifnot(beta<0)
    stopifnot(t>0)
    stopifnot(a<=b || b==Inf)
    stopifnot(a>=0)
    if(a==b) return(0)
    ##
    p1 <- if (b==Inf) Inf else (1.5*beta^2*t-log(b))/(beta*sqrt(t))
    p2 <- if (a==0) -Inf else if (a==Inf) Inf else (1.5*beta^2*t-log(a))/(beta*sqrt(t))
    ##
    temp <- pnorm(p1) - pnorm(p2)
    ##
    temp*exp(beta^2*t)
}

# Proportions -----------------------------------------------------------------------------------------------------

## see Appendix A2

zf_z <- function(t, z, beta) {
    x <- (0.5*beta^2*t - log(z))/(beta*sqrt(t))
    return (-dnorm(x)/(beta*sqrt(t)))
}

zh_z <- function(t, z, beta) {
    x <- (1.5*beta^2*t - log(z))/(beta*sqrt(t))
    p1 <- -dnorm(x)/(beta*sqrt(t))
    p2 <- pnorm(x)
    
    return (exp(beta^2*t)/z * (p1 - p2))
}

m <- function(t,z,beta){
    z*exp(beta^2*t)
}


# Find root -------------------------------------------------------------------------------------------------------

find_root_bi <- function(f,a,b,eps = sqrt(.Machine$double.eps)){
    ##
    stopifnot(a<b)
    stopifnot(eps>0)
    fa <- f(a)
    fb <- f(b)
    stopifnot(fa*fb<0)
    ##
    if ((b-a)<eps) return(0.5*(a+b)) # terminate condition
    ##
    c <- 0.5*(a + b)
    fc <- f(c)
    if (fc ==0) return(c) #impossible...
    ##
    if (fa*fc>0) {a <- c}
    if (fb*fc>0) {b <- c}
    ##
    find_root_bi(f,a,b,eps) # recursive step
}


# Determine conditions --------------------------------------------------------------------------------------------

Pr_for_U <- function(lambda, t, k, C, beta){
    stopifnot(lambda>0)
    stopifnot(t>0)
    stopifnot(k>C)
    stopifnot(beta<0)
    th <- (k-C)/lambda
    pnorm((log(th)+0.5*beta^2*t)/(-beta*sqrt(t)))
}
ES_for_U <- function(lambda, t, k, C, beta){
    stopifnot(lambda>0)
    stopifnot(t>0)
    stopifnot(k>C)
    stopifnot(beta<0)
    th <- (k-C)/lambda
    (C-k)*EZ0ab(th,Inf,t,beta)+lambda*EZ1ab(th,Inf,t,beta)
}
Q_ES_for_U <- function(lambda, t, k, C, beta){
    stopifnot(lambda>0)
    stopifnot(t>0)
    stopifnot(k>C)
    stopifnot(beta<0)
    th <- (k-C)/lambda
    (C-k)*EZ1ab(th,Inf,t,beta)+lambda*EZ2ab(th,Inf,t,beta)
}

