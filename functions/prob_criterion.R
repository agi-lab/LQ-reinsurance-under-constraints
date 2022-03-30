
# parameter -------------------------------------------------------------------------------------------------------


f_g1 <- function(
    g1,
    g2,
    T_, x, beta,k,
    C
){
    stopifnot(g1>0)
    stopifnot(g1<=g2)
    lambda <- (k-C)/g1
    ##
    E1 <- k*EZ1ab(0,g1,T_,beta) - lambda*EZ2ab(0,g1,T_,beta)
    E2 <- C*EZ1ab(g1,g2,T_,beta)
    E3 <- k*EZ1ab(g2,Inf,T_,beta) - lambda*EZ2ab(g2,Inf,T_,beta)
    ##
    E1+E2+E3-x
}


# path ------------------------------------------------------------------------------------------------------------

Xpi_P <- function(
    Z,t,
    lambda,T_,beta,a,b,k,
    C,
    g1,g2
){
    ## Z is Z_t
    stopifnot(Z>0)
    stopifnot(t>=0)
    ##
    lambda <- lambda*Z
    g1 <- g1/Z
    g2 <- g2/Z
    ##
    E1 <- k*EZ1ab(0,g1,T_-t,beta) - lambda*EZ2ab(0,g1,T_-t,beta)
    E2 <- C*EZ1ab(g1,g2,T_-t,beta)
    E3 <- k*EZ1ab(g2,Inf,T_-t,beta) - lambda*EZ2ab(g2,Inf,T_-t,beta)
    ##
    Xpi_tilda <- E1+E2+E3
    ## add back adjustment
    Xpi_tilda + (a-b)*t
}

pi_P <- function(Z, t,
                 lambda = NULL, T_, beta, sigma, k, C, 
                 g1, g2) {
    ## Z is Z_t
    ## lambda not used, inferred by g1 instead
    stopifnot(Z>0)
    stopifnot(t>=0)
    
    # c_P <- k - lambda*g2
    
    E1 <- k * (-zf_z(T_-t, g2/Z, beta) + zf_z(T_-t, g1/Z, beta))
    E2 <- (k-C) * ( g2/g1*zh_z(T_-t, g2/Z, beta) - zh_z(T_-t, g1/Z, beta)) + (k-C)/g1*m(T_-t,Z,beta)
    E3 <- C * (zf_z(T_-t, g2/Z, beta) - zf_z(T_-t, g1/Z, beta))
    
    return (1 + beta/sigma * (E1 + E2 + E3))
}
