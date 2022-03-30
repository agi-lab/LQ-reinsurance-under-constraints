source("functions/utils.R")

# Unconstrained
par_unconstrained <- function(a, b, sigma, T_, x, k0, C0) {
    beta = -1*b/sigma
    ## 
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    return (list(lambda_U = (k-x)*exp(-beta^2*T_)))
}

# Forced Constraint
par_constrained <- function(a, b, sigma, T_, x, k0, C0) {
    beta = -1*b/sigma
    ## 
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    # calculate lambda
    source("functions/forced_survival.R")
    ##
    f_lambdaC_loc <- function(lambda){
        f_lambdaC(lambda, T_=T_, beta = beta, x=x, k=k, C=C)
    }
    lambda_C <- find_root_bi(f_lambdaC_loc, 0.01, 100)
    
    return (list(lambda_C = lambda_C))
}

# Probability Constraint
par_prob_constraint <- function(a, b, sigma, T_, x, k0, C0, epsilon) {
    source("functions/prob_criterion.R")
    
    beta = -1*b/sigma
    ## 
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    
    Pr_U <- Pr_for_U(lambda_U, T_, k, C, beta)
    check_prob <- Pr_U > 1-epsilon
    stopifnot(!check_prob)
    
    # get parameters
    g2 <- exp(-0.5*beta^2*T_-beta*sqrt(T_)*qnorm(1-epsilon))
    f_g1_loc <- function(g1){
        f_g1(g1,g2=g2,T_=T_, x=x, beta=beta,k=k,C=C)
    }
    ##
    g1 <- find_root_bi(f_g1_loc, 0.01, g2)
    ## 
    lambda_P <- (k-C)/g1
    c_P <- k-lambda_P*g2
    
    return (list(lambda_P = lambda_P, c_P = c_P, g1 = g1, g2 = g2))
}
    
# Expected Shortfall (P)
par_ES_constraint_P <- function(a, b, sigma, T_, x, k0, C0, nu) {
    source("functions/expected_shortfall.R")
    
    beta = -1*b/sigma
    ## 
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    
    ES_U <- ES_for_U(lambda_U, T_, k, C, beta)
    check_ES <- ES_U < nu
    stopifnot(!check_ES)
    
    f_h0_loc <- function(h){
        f_h0(h, k=k, C=C, nu=nu, beta = beta, t=T_)
    }
    h0 <- find_root_bi(f_h0_loc,0.01,100)
    
    f_h_loc <- function(h){
        f_h(h,h0=h0, x=x, k=k,C=C, nu=nu, beta = beta, t=T_)
    }
    ## parameters
    h2 <- find_root_bi(f_h_loc, h0*1.01,100)
    lambda_E <- f_lambdaP(h2, nu, beta=beta, t=T_)
    h1 <- (k-C)/lambda_E
    gamma_E <- (h2-h1)*lambda_E
    
    return (list(lambda_E = lambda_E, gamma_E = gamma_E, h1 = h1, h2 = h2))
}

# Expected Shortfall (Q)
    
par_ES_constraint_Q <- function(a, b, sigma, T_, x, k0, C0, nu) {
    source("functions/Q_shortfall.R")
    
    beta = -1*b/sigma
    ## 
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    
    Q_ES_U <- Q_ES_for_U(lambda_U, T_, k, C, beta)
    check_Q_ES <- Q_ES_U < nu
    check_Q_ES
    
    f_delta_loc <- function(delta){
        f_delta(delta, k=k,C=C,nu=nu,t=T_,beta=beta)
    }
    delta_Q <- find_root_bi(f_delta_loc,0,10)
    f_lambdaQ_loc <- function(lambda){
        f_lambdaQ(lambda, x=x, k=k,C=C,delta=delta_Q,t=T_,beta=beta)
    }
    lambda_Q <- find_root_bi(f_lambdaQ_loc,delta_Q,100)
    
    return (list(lambda_Q = lambda_Q, delta_Q = delta_Q))
}
