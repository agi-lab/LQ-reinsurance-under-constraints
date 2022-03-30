
# parameters ------------------------------------------------------------------------------------------------------


f_lambda_h <- function(h, ...){
    ## `...` for t, beta
    stopifnot(h>=0)
    ##
    EZ1ab(h,Inf, ...) - h*EZ0ab(h,Inf, ...)
}
f_lambdaP <- function(h, nu, ...){
    nu / f_lambda_h(h, ...)
}

f_h0 <- function(h, k, C, nu, ...){
    lambda <- f_lambdaP(h, nu, ...)
    h1 <- (k-C)/lambda
    ##
    h-h1
}

f_h <- function(
    h,
    h0,
    x,
    nu,
    k,C,
    ...
){
    stopifnot(h>h0)
    lambda <- f_lambdaP(h, nu, ...)
    h1 <- (k-C)/lambda
    h2 <- h
    gamma <- (h2-h1)*lambda
    ##
    E1 <- k*EZ1ab(0,h1,...) - lambda*EZ2ab(0,h1,...)
    E2 <- C*EZ1ab(h1,h2,...)
    E3 <- (k+gamma)*EZ1ab(h2,Inf,...) - lambda*EZ2ab(h2,Inf,...)
    ##
    E1+E2+E3-x
}


# path ------------------------------------------------------------------------------------------------------------

Xpi_E <- function(
    Z,t,
    lambda,T_,beta,a,b,k,
    C,
    h1,h2
){
    ## Z is Z_t
    stopifnot(Z>0)
    stopifnot(t>=0)
    gamma = (h2-h1)*lambda
    ##
    lambda <- lambda*Z
    h1 <- h1/Z
    h2 <- h2/Z
    ##
    E1 <- k*EZ1ab(0,h1,T_-t,beta) - lambda*EZ2ab(0,h1,T_-t,beta)
    E2 <- C*EZ1ab(h1,h2,T_-t,beta)
    E3 <- (k+gamma)*EZ1ab(h2,Inf,T_-t,beta) - lambda*EZ2ab(h2,Inf,T_-t,beta)
    ##
    E1+E2+E3
    ##
    Xpi_tilda <- E1+E2+E3
    ## add back adjustment
    Xpi_tilda + (a-b)*t
}

pi_E <- function(Z, t,
                 lambda=NULL, T_, beta, sigma, k, C,
                 h1, h2) {
    ## Z is Z_t
    ## lambda not used, inferred by h1 instead
    stopifnot(Z>0)
    stopifnot(t>=0)
    
    gamma <- (h2 - h1)*(k - C)/h1
    
    E1 <- k * (-zf_z(T_-t, h2/Z, beta) + zf_z(T_-t, h1/Z, beta))
    E2 <- (k-C) * (h2/h1*zh_z(T_-t, h2/Z, beta) - zh_z(T_-t, h1/Z, beta))  + (k-C)/h1*m(T_-t,Z,beta)
    E3 <- C * (zf_z(T_-t, h2/Z, beta) - zf_z(T_-t, h1/Z, beta))
    E4 <- gamma * -zf_z(T_-t, h2/Z, beta)
        
    return (1 + beta/sigma * (E1 + E2 + E3 + E4))
}

