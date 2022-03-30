
# parameters ------------------------------------------------------------------------------------------------------

f_delta <- function(delta, k, C, nu, ...){
    ## `..` for t, beta
    stopifnot(delta>=0)
    th <- if(delta==0) Inf else (k-C)/delta
    ##
    delta*EZ2ab(th,Inf,...) + (C-k)*EZ1ab(th,Inf,...)-nu
}
f_lambdaQ <- function(
    lambda, x,
    k, C, delta,
    ...
){
    stopifnot(lambda>=delta)
    ##
    a1 <- (k-C)/lambda
    b1 <- (k-C)/delta
    ##
    E1 <- k*EZ1ab(0,a1,...)-lambda*EZ2ab(0,a1,...)
    E2 <- C*EZ1ab(a1,b1,...)
    E3 <- k*EZ1ab(b1,Inf,...)-delta*EZ2ab(b1,Inf,...)
    ##
    E1+E2+E3-x
}

# path ------------------------------------------------------------------------------------------------------------

Xpi_Q <- function(
    Z,t,
    T_,beta,a,b,k,
    C,
    lambda, delta
){
    ## Z is Z_t
    stopifnot(Z>0)
    stopifnot(t>=0)
    ##
    lambda <- lambda*Z
    delta <- delta*Z
    a1 <- (k-C)/lambda
    b1 <- (k-C)/delta
    ##
    E1 <- k*EZ1ab(0,a1,T_-t,beta)-lambda*EZ2ab(0,a1,T_-t,beta)
    E2 <- C*EZ1ab(a1,b1,T_-t,beta)
    E3 <- k*EZ1ab(b1,Inf,T_-t,beta)-delta*EZ2ab(b1,Inf,T_-t,beta)
    ##
    Xpi_tilda <- E1+E2+E3
    ## add back adjustment
    Xpi_tilda + (a-b)*t
}

pi_Q <- function(Z, t,
                 lambda, T_, beta, sigma, k, C,
                 delta) {
    ## Z is Z_t
    stopifnot(Z>0)
    stopifnot(t>=0)
    
    a_lambda <- (k-C)/lambda
    a_delta <- (k-C)/delta
    
    E1 <- k * (zf_z(T_-t, a_lambda/Z, beta) - zf_z(T_-t, a_delta/Z, beta))
    E2 <- C * (zf_z(T_-t, a_delta/Z, beta) - zf_z(T_-t, a_lambda/Z, beta))
    E3 <- -(k-C) * zh_z(T_-t, a_lambda/Z, beta)
    E4 <- (k-C) * zh_z(T_-t, a_delta/Z, beta) + delta*m(T_-t,Z,beta)
        
    return (1 + beta/sigma * (E1 + E2 + E3 + E4))
}
