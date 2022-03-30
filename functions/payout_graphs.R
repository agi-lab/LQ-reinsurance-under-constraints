library(ggplot2)
source("functions/get_parameters.R")
source("functions/utils.R")

# No Constraint
payout_unconstrained <- function(a, b, sigma, T_, x, k0, C0, lwd=1.2) {

    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    
    Z_U <- (k - C)/lambda_U
    
    start_point <- c(-2*Z_U, k - lambda_U*2*Z_U)
    end_point <- c(0, k)
    mid_point = 0.5*(start_point + end_point)

    gplot <- ggplot() +
        geom_segment(aes(x = start_point[1], xend = end_point[1], 
                         y = start_point[2], yend = end_point[2]),size=lwd) +
        geom_hline(aes(yintercept=C),linetype='dashed') +
        geom_vline(aes(xintercept=-Z_U),linetype='dashed') +
        annotate(geom = "text",
                 x = 0, y = k,
                 label = "k",
                 hjust = 0, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = mid_point[1]*0.5, y = C,
                 label = "italic(C)",
                 hjust = 0, vjust = 1.5, parse=TRUE) +
        annotate(geom = "text",
                 x = -Z_U, y = 0.5*(mid_point[2] + end_point[2]),
                 label = "italic((k - C)/lambda[U])",
                 hjust = 1.2, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = 0.5*(mid_point[1]+start_point[1]), 
                 y = 0.5*(mid_point[2]+start_point[2]),
                 label = "Unconstrained Payoff",
                 hjust = 0, vjust = 1.2) +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks = element_blank())+
        xlab("") + ylab("")
    
    return (gplot)
}

# Forced Constraint
payout_forced_constraint <- function(a, b, sigma, T_, x, k0, C0, lwd=1.2) {

    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_

    # Find lambda_C
    lambda_C <- par_constrained(a, b, sigma, T_, x, k0, C0)$lambda_C
    Z_C <- (k - C)/lambda_C
    
    start_point_C <- c(-Z_C, k - lambda_C*Z_C)
    end_point_C <- c(0, k)
    mid_point_C <- c(-Z_C, C)
        
    #Find lambda_U
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    Z_U <- (k - C)/lambda_U
    start_point_U <- c(-Z_U-Z_C, k - lambda_U*(Z_U+Z_C))
    end_point_U <- c(0, k)
    
    gplot <- ggplot() +
        # Above C
        geom_segment(aes(x = start_point_C[1], xend = end_point_C[1], 
                         y = start_point_C[2], yend = end_point_C[2]), size=lwd) +
        # Max C
        geom_segment(aes(x = -Z_U-Z_C, xend = -Z_C, 
                         y = C, yend = C),size=lwd) +
        # Baseline
        geom_segment(aes(x = start_point_U[1], xend = end_point_U[1],
                         y = start_point_U[2], yend = end_point_U[2]), 
                     size=1, linetype = 'dashed', col = 'darkgrey') +
        geom_hline(aes(yintercept=C),linetype='dashed') +
        geom_vline(aes(xintercept=-Z_C),linetype='dashed') +
        annotate(geom = "text",
                 x = 0, y = k,
                 label = "k",
                 hjust = 0, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = 0.5*(start_point_C[1]+end_point_C[1]), y = C,
                 label = "italic(C)",
                 hjust = 0, vjust = 1.2, parse=TRUE) +
        annotate(geom = "text",
                 x = -Z_C, y = 0.5*(C + start_point_U[2]),
                 label = "italic((k - C)/lambda[C])",
                 hjust = 1.2, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = 0.5*(start_point_U[1]+end_point_U[1]),
                 y = 0.5*(start_point_U[2]+end_point_U[2]),
                 label = "Unconstrained Payoff",
                 hjust = 0, vjust = 1.2) +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks = element_blank())+
        xlab("") + ylab("")
    
    return (gplot)
}

# Probability Constraint
payout_prob_constraint <- function(a, b, sigma, T_, x, k0, C0, epsilon, lwd=1.2) {
    
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    prob_par <- par_prob_constraint(a, b, sigma, T_, x, k0, C0, epsilon)
    lambda_P <- prob_par$lambda_P
    c_P <- prob_par$c_P
    
    Z_c <- (k - c_P)/lambda_P
    Z_C <- (k - C)/lambda_P
    
    #Find lambda_U
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    Z_U <- (k - C)/lambda_U
    start_point_U <- c(-(Z_c+Z_C), k - lambda_U*(Z_c+Z_C))
    end_point_U <- c(0, k)
    
    gplot <- ggplot() +
        # Above C
        geom_segment(aes(x = -Z_C, xend = 0, 
                         y = C, yend = k),size=lwd) +
        # Below c_P
        geom_segment(aes(x = -(Z_c+Z_C), xend = -Z_c, 
                         y = k - lambda_P*(Z_c+Z_C), yend = c_P),size=lwd) +
        # Between c and C
        geom_segment(aes(x = -Z_c, xend = -Z_C, 
                         y = C, yend = C),size=lwd) +
        # Baseline
        geom_segment(aes(x = start_point_U[1], xend = end_point_U[1],
                         y = start_point_U[2], yend = end_point_U[2]),
                     size=1, linetype = 'dashed', col = 'darkgrey') +
        geom_hline(aes(yintercept=C),linetype='dashed') +
        geom_hline(aes(yintercept=c_P),linetype='dashed') +
        geom_vline(aes(xintercept=-Z_C),linetype='dashed') +
        geom_vline(aes(xintercept=-Z_c),linetype='dashed') +
        annotate(geom = "text",
                 x = 0, y = k,
                 label = "k",
                 hjust = 0, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = -Z_C*0.5, y = C,
                 label = "italic(C)",
                 hjust = 0, vjust = 1.5, parse=TRUE) +
        annotate(geom = "text",
                 x = -Z_C*0.5, y = c_P,
                 label = "italic(c[P])",
                 hjust = 0, vjust = 1.5, parse=TRUE) +
        annotate(geom = "text",
                 x = -Z_C, y = 0.5*(C + start_point_U[2]),
                 label = "italic(g[1] : (k - C)/lambda[P])",
                 hjust = 1.2, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = -Z_c, y = 0.5*(C + start_point_U[2]),
                 label = "italic(g[2] : (k - c[P])/lambda[P])",
                 hjust = 1.2, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = 0.5*(start_point_U[1]+end_point_U[1]),
                 y = 0.5*(start_point_U[2]+end_point_U[2]),
                 label = "Unconstrained Payoff",
                 hjust = 0, vjust = 1.2) +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks = element_blank())+
        xlab("") + ylab("")
    
    return (gplot)
}

# Expected Shortfall Constraint (P)
payout_ES_constraint_P <- function(a, b, sigma, T_, x, k0, C0, nu, lwd=1.2) {
    
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    ES_P_par <- par_ES_constraint_P(a, b, sigma, T_, x, k0, C0, nu)
    lambda_E <- ES_P_par$lambda_E
    gamma_E <- ES_P_par$gamma_E
    
    Z_C <- (k-C)/lambda_E
    Z_gamma <- (k - (C - gamma_E))/lambda_E
    
    #Find lambda_U
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    Z_U <- (k - C)/lambda_U
    start_point_U <- c(-Z_U-Z_gamma, k - lambda_U*(Z_U+Z_gamma))
    end_point_U <- c(0, k)
    
    gplot <- ggplot() +
        # Above C
        geom_segment(aes(x = -Z_C, xend = 0, 
                         y = C, yend = k),size=lwd) +
        # Below C - gamma_E
        geom_segment(aes(x = -(Z_gamma+Z_U), xend = -Z_gamma,
                         y = k - lambda_E*(Z_gamma+Z_U) + gamma_E, yend = k - lambda_E*Z_gamma + gamma_E),size=lwd) +
        # Between C - gamma_E and C
        geom_segment(aes(x = -Z_gamma, xend = -Z_C,
                         y = C, yend = C),size=lwd) +
        # # Baseline
        geom_segment(aes(x = start_point_U[1], xend = end_point_U[1],
                         y = start_point_U[2], yend = end_point_U[2]),
                     size=1, linetype = 'dashed', col = 'darkgrey') +
        geom_hline(aes(yintercept=C),linetype='dashed') +
        geom_hline(aes(yintercept=C - gamma_E),linetype='dashed') +
        geom_vline(aes(xintercept=-Z_C),linetype='dashed') +
        geom_vline(aes(xintercept=-Z_gamma),linetype='dashed') +
        annotate(geom = "text",
                 x = 0, y = k,
                 label = "k",
                 hjust = 0, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = -0.5*Z_C, y = C,
                 label = "italic(C)",
                 hjust = 0, vjust = 1.5, parse=TRUE) +
        annotate(geom = "text",
                 x = -0.5*Z_C, y = C - gamma_E,
                 label = "italic(C - gamma[S])",
                 hjust = 0, vjust = 1.5, parse=TRUE) +
        annotate(geom = "text",
                 x = -Z_C, y = 0.5*(k+C),
                 label = "italic(h[1] : (k-C)/lambda[S])",
                 hjust = 1.1, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = -Z_gamma, y = 0.5*(k+C),
                 label = "italic(h[2] : (k + gamma[S]-C)/lambda[S])",
                 hjust = 1.1, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = 0.5*(start_point_U[1]+end_point_U[1]),
                 y = 0.5*(start_point_U[2]+end_point_U[2]),
                 label = "Unconstrained Payoff",
                 hjust = 0, vjust = 1.2) +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks = element_blank())+
        xlab("") + ylab("")
    
    return (gplot)
}

# Expected Shortfall Constraint (Q)

payout_ES_constraint_Q <- function(a, b, sigma, T_, x, k0, C0, nu, lwd=1.2) {
    
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    ES_Q_par <- par_ES_constraint_Q(a, b, sigma, T_, x, k0, C0, nu)
    lambda_Q <- ES_Q_par$lambda_Q
    delta_Q <- ES_Q_par$delta_Q
    
    ES_Q_interval <- (k-C)/c(lambda_Q, delta_Q)
    
    #Find lambda_U
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    Z_U <- (k - C)/lambda_U
    start_point_U <- c(-(ES_Q_interval[1] + ES_Q_interval[2]), k - lambda_U*(ES_Q_interval[1] + ES_Q_interval[2]))
    end_point_U <- c(0, k)
    
    gplot <- ggplot() +
        # Above upper
        geom_segment(aes(x = -ES_Q_interval[1], xend = 0, 
                         y = k - lambda_Q*ES_Q_interval[1], yend = k),size=lwd) +
        # Below lower
        geom_segment(aes(x = -(ES_Q_interval[2] + ES_Q_interval[1]), xend = -ES_Q_interval[2],
                         y = k - delta_Q*(ES_Q_interval[2] + ES_Q_interval[1]), yend = k - delta_Q*(ES_Q_interval[2])),size=lwd) +
        # Between upper and lower
        geom_segment(aes(x = -ES_Q_interval[2], xend = -ES_Q_interval[1],
                         y = C, yend = C),size=lwd) +
        # # Baseline
        geom_segment(aes(x = start_point_U[1], xend = end_point_U[1],
                         y = start_point_U[2], yend = end_point_U[2]),
                     size=1, linetype = 'dashed', col = 'darkgrey') +
        # guiding line
        geom_segment(aes(x = -ES_Q_interval[2], xend = end_point_U[1],
                         y = C, yend = k),
                     size=1, linetype = 'dashed', col = 'lightgrey') +
        geom_hline(aes(yintercept=C),linetype='dashed') +
        geom_vline(aes(xintercept=-ES_Q_interval[1]),linetype='dashed') +
        geom_vline(aes(xintercept=-ES_Q_interval[2]),linetype='dashed') +
        annotate(geom = "text",
                 x = 0, y = k,
                 label = "k",
                 hjust = 0, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = -0.5*ES_Q_interval[1], y = C,
                 label = "italic(C)",
                 hjust = 0, vjust = 1.5, parse=TRUE) +
        annotate(geom = "text",
                 x = -ES_Q_interval[1], y = 0.5*(k+C),
                 label = "italic((k-C)/lambda[Q])",
                 hjust = 1.1, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = -ES_Q_interval[2], y = 0.5*(k+C),
                 label = "italic((k - C)/delta[Q])",
                 hjust = -0.1, vjust = 0, parse=TRUE) +
        annotate(geom = "text",
                 x = 0.5*(start_point_U[1]+end_point_U[1]),
                 y = 0.5*(start_point_U[2]+end_point_U[2]),
                 label = "Unconstrained Payoff",
                 hjust = 0, vjust = 1.2) +
        coord_cartesian(ylim = c(C-k, k)) +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks = element_blank())+
        xlab("") + ylab("")
    return (gplot)
}

payout_graphs <- function(a, b, sigma, T_, x, k0, C0, epsilon = 0, nu = 0,
                          UC=F, FC=F, PC=F, ESC_P=F, ESC_Q=F, lwd = 1.2) {
    
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    gplot <- ggplot()
    x.lim <- c(0, 0)
    y.lim <- c(0, k)
    
    scale_breaks = c()
    scale_values = c()
    
    if (FC) {
        # Handling Legends
        FC_name <- "Strict Constraint"
        FC_col <- "blue"
        
        scale_breaks <- c(scale_breaks, FC_name)
        scale_values <- c(scale_values, "Strict Constraint" = FC_col)
        
        # Calculations
        lambda_C <- par_constrained(a, b, sigma, T_, x, k0, C0)$lambda_C
        Z_C <- (k - C)/lambda_C
        
        x.lim <- c(min(x.lim[1], -2.1*Z_C), max(x.lim[2], 0.1))
        y.lim <- c(min(y.lim[1], -1.1*k), max(y.lim[2], 1.1*k))
        
        gplot <- gplot +
            # Above Z_C
            geom_segment(aes(x = -Z_C, xend = 0, 
                             y = k - lambda_C*Z_C, yend = k, color = FC_name), size=lwd) +
            # Below Z_C
            geom_segment(aes(x = -Z_C*1000, xend = -Z_C, 
                             y = C, yend = C, color = FC_name), size=lwd) +
            geom_vline(aes(xintercept=-Z_C),linetype='dashed', col = "dodgerblue1") +
            annotate(geom = "text",
                     x = -Z_C, y = mean(y.lim),
                     label = "italic((k - C)/lambda[C])",
                     hjust = -0.1, vjust = 0, parse=TRUE, col = "dodgerblue1")
        
    }
    if (PC) {
        stopifnot(epsilon > 0)
        # Handling legend
        PC_name <- "VaR Constraint"
        PC_col <- "green4"
        
        scale_breaks <- c(scale_breaks, PC_name)
        scale_values <- c(scale_values, "VaR Constraint" = PC_col)
        
        # Perform calculations
        prob_par <- par_prob_constraint(a, b, sigma, T_, x, k0, C0, epsilon)
        lambda_P <- prob_par$lambda_P
        c_P <- prob_par$c_P
        
        Z_CP <- (k - C)/lambda_P
        Z_cP <- (k - c_P)/lambda_P
        
        x.lim <- c(min(x.lim[1], -3/2*Z_cP), max(x.lim[2], 0.1))
        y.lim <- c(min(y.lim[1], -1.1*k), max(y.lim[2], 1.1*k))
        
        gplot <- gplot +
            # Above C
            geom_segment(aes(x = -Z_CP, xend = 0, 
                             y = C, yend = k, color = PC_name), size=lwd) +
            # Below c_P
            geom_segment(aes(x = -(1000*Z_cP), xend = -Z_cP, 
                             y = k - lambda_P*(1000*Z_cP), yend = c_P, color = PC_name), size=lwd) +
            # Between c and C
            geom_segment(aes(x = -Z_cP, xend = -Z_CP, 
                             y = C, yend = C, color = PC_name), size=lwd) +
            geom_hline(aes(yintercept=c_P),linetype='dashed', col = "seagreen4") +
            geom_vline(aes(xintercept=-Z_cP),linetype='dashed', col = "seagreen4") +
            annotate(geom = "text",
                     x = -0.5*(Z_CP + Z_cP), y = c_P,
                     label = "italic(c[P])",
                     hjust = 0, vjust = 1.5, parse=TRUE, col = "seagreen4") +
            annotate(geom = "text",
                     x = -Z_cP, y = k,
                     label = "italic((k - c[P])/lambda[P])",
                     hjust = 1.2, vjust = 0, parse=TRUE, col = "seagreen4")
        
    }
    if (ESC_P) {
        stopifnot(nu > 0)
        # Handling legend
        ESC_P_name <- "P-ES Constraint"
        ESC_P_col <- "red3"
        
        scale_breaks <- c(scale_breaks, ESC_P_name)
        scale_values <- c(scale_values, "P-ES Constraint" = ESC_P_col)
        
        # Perform calculations
        ES_P_par <- par_ES_constraint_P(a, b, sigma, T_, x, k0, C0, nu)
        lambda_E <- ES_P_par$lambda_E
        gamma_E <- ES_P_par$gamma_E
        
        Z_CE <- (k - C)/lambda_E
        Z_CgE <- (k - (C - gamma_E))/lambda_E
        
        x.lim <- c(min(x.lim[1], -3/2*Z_CgE), max(x.lim[2], 0.1))
        y.lim <- c(min(y.lim[1], -1.1*k), max(y.lim[2], 1.1*k))
        
        gplot <- gplot +
            # Above C
            geom_segment(aes(x = -Z_CE, xend = 0, 
                             y = C, yend = k, color = ESC_P_name), size=lwd) +
            # Below C - gamma_E
            geom_segment(aes(x = -(1000*Z_CgE), xend = -Z_CgE,
                             y = k - lambda_E*(1000*Z_CgE) + gamma_E, yend = k - lambda_E*Z_CgE + gamma_E,
                             color = ESC_P_name), size=lwd) +
            # Between C - gamma_E and C
            geom_segment(aes(x = -Z_CgE, xend = -Z_CE,
                             y = C, yend = C, color = ESC_P_name), size=lwd) +
            geom_vline(aes(xintercept=-Z_CE),linetype='dashed', col = "tomato2") +
            geom_vline(aes(xintercept=-Z_CgE),linetype='dashed', col = "tomato2") +
            annotate(geom = "text",
                     x = -Z_CE, y = 0.5*(k+C),
                     label = "italic((k-C)/lambda[S])",
                     hjust = 1.1, vjust = 0, parse=TRUE, col = "tomato2") +
            annotate(geom = "text",
                     x = -Z_CgE, y = 0.5*(k+C),
                     label = "italic((k + gamma[S] - C)/lambda[S])",
                     hjust = 1.1, vjust = 0, parse=TRUE, col = "tomato2")
        
    }
    if (ESC_Q) {
        stopifnot(nu > 0)
        # Handling legend
        ESC_Q_name <- "Q-ES Constraint"
        ESC_Q_col <- "orchid4"
        
        scale_breaks <- c(scale_breaks, ESC_Q_name)
        scale_values <- c(scale_values, "Q-ES Constraint" = ESC_Q_col)
        
        # Perform calculations
        ES_Q_par <- par_ES_constraint_Q(a, b, sigma, T_, x, k0, C0, nu)
        lambda_Q <- ES_Q_par$lambda_Q
        delta_Q <- ES_Q_par$delta_Q
        
        Z_lQ <- (k - C)/lambda_Q
        Z_ld <- (k - C)/delta_Q
        
        x.lim <- c(min(x.lim[1], -3/2*Z_ld), max(x.lim[2], 0.1))
        y.lim <- c(min(y.lim[1], -1.1*k), max(y.lim[2], 1.1*k))
        
        gplot <- gplot +
            # Above upper
            geom_segment(aes(x = -Z_lQ, xend = 0, 
                             y = C, yend = k, color = ESC_Q_name), size=lwd) +
            # Below lower
            geom_segment(aes(x = -(1000*Z_ld), xend = -Z_ld,
                             y = k - delta_Q*(1000*Z_ld), yend = k - delta_Q*(Z_ld), color = ESC_Q_name), size=lwd) +
            # Between upper and lower
            geom_segment(aes(x = -Z_ld, xend = -Z_lQ,
                             y = C, yend = C, color = ESC_Q_name), size=lwd) +
            geom_vline(aes(xintercept=-Z_lQ),linetype='dashed', col = "mediumpurple2") +
            geom_vline(aes(xintercept=-Z_ld),linetype='dashed', col = "mediumpurple2") +
            annotate(geom = "text",
                     x = -Z_lQ, y = k,
                     label = "italic((k-C)/lambda[Q])",
                     hjust = 1.1, vjust = 0, parse=TRUE, col = "mediumpurple2") +
            annotate(geom = "text",
                     x = -Z_ld, y = k,
                     label = "italic((k - C)/delta[Q])",
                     hjust = 1.1, vjust = 0, parse=TRUE, col = "mediumpurple2")
        
    }
    if (UC) {
        # Handling legend
        UC_name <- "Unconstrained"
        UC_col <- "darkgrey"
        
        scale_breaks <- c(scale_breaks, UC_name)
        scale_values <- c(scale_values, "Unconstrained" = UC_col)
        
        # Perform calculations
        #Find lambda_U
        lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
        Z_U <- (k - C)/lambda_U
        
        x.lim <- c(min(x.lim[1], -2.1*Z_U), max(x.lim[2], 0.1))
        y.lim <- c(min(y.lim[1], -1.1*k), max(y.lim[2], 1.1*k))
        
        gplot <- gplot + 
            geom_segment(aes(x = -1000*Z_U, xend = 0, 
                             y = k - lambda_U*1000*Z_U, yend = k, color = UC_name),
                         linetype='dashed', size=lwd)
    }
    
    gplot <- gplot +
        geom_hline(aes(yintercept=C),linetype='dashed') +
        annotate(geom = "text",
                 x = 0, y = C,
                 label = "italic(C)",
                 hjust = 0, vjust = 1.2, parse=TRUE) +
        annotate(geom = "text",
                 x = 0, y = k,
                 label = "k",
                 hjust = -1, vjust = 0, parse=TRUE) +
        coord_cartesian(xlim = x.lim, ylim = y.lim) +
        scale_colour_manual(name = "Constraints",
                            breaks = scale_breaks,
                            values = scale_values) +
        # theme(panel.grid.major = element_blank(), 
        #       panel.grid.minor = element_blank(),
        #       panel.background = element_blank(),
        #       axis.text.x=element_blank(),
        #       axis.text.y=element_blank(),
        #       axis.ticks = element_blank()) +
        xlab("") + ylab("")

    return (gplot)
}

######################################################
# payout_unconstrained(0.2, 0.5, 1.2, 5, 2, 5, 0)
# payout_forced_constraint(0.2, 0.5, 1.2, 5, 2, 5, 0)
# payout_prob_constraint(0.2, 0.5, 1.2, 5, 2, 5, 0, 0.01)
# payout_ES_constraint_P(0.2, 0.5, 1.2, 5, 2, 5, 0, 0.1)
# payout_ES_constraint_Q(0.2, 0.5, 1.2, 5, 2, 5, 0, 0.1)
# 
# payout_graphs(0.2, 0.5, 1.2, 5, 2, 5, 0,
#               epsilon = 0.01, nu = 0.1,
#               UC = T, FC = T, PC = T, ESC_P = T, ESC_Q = T)
# 
# payout_graphs(0.2, 0.5, 1.2, 5, 2, 5, 0,
#               epsilon = 0.01, nu = 0.1,
#               FC = T, ESC_Q=T)
# payout_graphs(0.2, 0.5, 1.2, 5, 2, 5, 0,
#               epsilon = 0.01, nu = 0.1,
#               PC = T, ESC_P=T)
# 
# 
# ######################################################
# # TO PDF
# 
# pdf("figures/unconstrained.pdf", height = 7)
# payout_unconstrained(0.2, 0.5, 1.2, 5, 2, 5, 0)
# dev.off()
# 
# pdf("figures/strict constraint.pdf", height = 7)
# payout_forced_constraint(0.2, 0.5, 1.2, 5, 2, 5, 0)
# dev.off()
# 
# pdf("figures/var constraint.pdf", height = 7)
# payout_prob_constraint(0.2, 0.5, 1.2, 5, 2, 5, 0, 0.01)
# dev.off()
# 
# pdf("figures/ES constraint in P.pdf", height = 7)
# payout_ES_constraint_P(0.2, 0.5, 1.2, 5, 2, 5, 0, 0.1)
# dev.off()
# 
# pdf("figures/ES constraint in Q.pdf", height = 7)
# payout_ES_constraint_Q(0.2, 0.5, 1.2, 5, 2, 5, 0, 0.1)
# dev.off()
# 
# pdf("figures/Payout SC and ESQ.pdf", height = 7)
# payout_graphs(0.2, 0.5, 1.2, 5, 2, 5, 0,
#               epsilon = 0.01, nu = 0.1,
#               FC = T, ESC_Q=T)
# dev.off()
# 
# pdf("figures/Payout VaR and ESP.pdf", height = 7)
# payout_graphs(0.2, 0.5, 1.2, 5, 2, 5, 0,
#               epsilon = 0.01, nu = 0.1,
#               PC = T, ESC_P=T)
# dev.off()