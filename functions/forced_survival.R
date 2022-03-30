
# parameter -------------------------------------------------------------------------------------------------------

f_lambdaC <- function(lambda, T_, x, beta, k, C){
    stopifnot(lambda>0)
    ##
    th <- (k-C)/lambda
    E1 <- k*EZ1ab(0,th,T_,beta) - lambda*EZ2ab(0,th,T_,beta)
    E2 <- C*EZ1ab(th,Inf,T_,beta)
    ##
    E1+E2-x 
}


# path ------------------------------------------------------------------------------------------------------------

Xpi_C <- function(
    Z,t,
    lambda,T_,beta,a,b,k,
    C 
){
    ## Z is Z_t
    stopifnot(Z>0)
    stopifnot(t>=0)
    ##
    lambda <- lambda*Z
    ##
    th <- (k-C)/lambda
    E1 <- k*EZ1ab(0,th,T_-t,beta) - lambda*EZ2ab(0,th,T_-t,beta)
    E2 <- C*EZ1ab(th,Inf,T_-t,beta)
    ##
    Xpi_tilda <- E1+E2
    ## add back adjustment
    Xpi_tilda + (a-b)*t
}

pi_C <- function(Z, t,
                 lambda, T_, beta, sigma, k, C) {
    ## Z is Z_t
    stopifnot(Z>0)
    stopifnot(t>=0)
    
    eta <- (k-C)/lambda
    
    E1 <- k * zf_z(T_-t, eta/Z, beta)
    E2 <- -(k-C) * zh_z(T_-t, eta/Z, beta)
    E3 <- -C *zf_z(T_-t, eta/Z, beta)
    
    return (1 + beta/sigma * (E1 + E2 + E3))
}
