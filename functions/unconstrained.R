
# path ------------------------------------------------------------------------------------------------------------

Xpi_U <- function(
    Z,t,
    lambda,T_,beta,a,b,k
){
    ## Z is Z_t
    stopifnot(Z>0)
    stopifnot(t>=0)
    lambda <- lambda*Z
    Xpi_tilda <- k - lambda*EZ2ab(0,Inf,T_-t,beta)
    ## add back adjustment
    Xpi_tilda + (a-b)*t
}

pi_U <- function(Z, t, 
                 lambda, T_, beta, sigma) {
    stopifnot(Z>0)
    stopifnot(t>=0)
    
    1 + beta/sigma*(lambda*m(T_-t,Z,beta))
}
