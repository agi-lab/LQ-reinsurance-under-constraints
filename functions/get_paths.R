library(withr)
library(purrr)
library(tidyr)
library(ggplot2)
library(cowplot)
source("functions/utils.R")
source("functions/get_parameters.R")
source("functions/unconstrained.R")
source("functions/forced_survival.R")
source("functions/prob_criterion.R")
source("functions/expected_shortfall.R")
source("functions/Q_shortfall.R")

# X path ------------------------------------------------------------------------------------------------------------

calc_paths_X <- function(seed, Nseq, 
                       a, b, sigma, T_, x, k0, C0, epsilon, nu,
                       tilde = T) {
    beta = -1*b/sigma
    ## 
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    ##
    
    tseq <- seq(0,T_,length.out = Nseq+1)[-1]
    tseq[length(tseq)] <- tseq[length(tseq)] - 0.00001
    dt <- tseq[2]-tseq[1]
    
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    Xp_U <- function(Z, t){
        Xpi_U(Z,t,lambda = lambda_U, T_=T_, beta = beta, a=a, b=b, k=k)
    }
    
    lambda_C <- par_constrained(a, b, sigma, T_, x, k0, C0)$lambda_C
    Xp_C <- function(Z, t){
        Xpi_C(Z,t,lambda = lambda_C, T_=T_, beta = beta, a=a, b=b, k=k, C=C)
    }
    
    par_P <- par_prob_constraint(a, b, sigma, T_, x, k0, C0, epsilon)
    lambda_P <- par_P$lambda_P
    c_P <- par_P$c_P
    g1 <- par_P$g1
    g2 <- par_P$g2
    Xp_P <- function(Z, t){
        Xpi_P(Z,t,lambda = lambda_P, T_=T_, beta = beta, a=a, b=b, k=k, C=C, g1=g1, g2=g2)
    }
    
    par_E <- par_ES_constraint_P(a, b, sigma, T_, x, k0, C0, nu)
    lambda_E <- par_E$lambda_E
    gamma_E <- par_E$gamma_E
    h1 <- par_E$h1
    h2 <- par_E$h2
    Xp_E <- function(Z, t){
        Xpi_E(Z,t,lambda = lambda_E, T_=T_, beta = beta, a=a, b=b, k=k, C=C, h1=h1, h2=h2)
    }
    
    par_Q <- par_ES_constraint_Q(a, b, sigma, T_, x, k0, C0, nu)
    lambda_Q <- par_Q$lambda_Q
    delta_Q <- par_Q$delta_Q
    Xp_Q <- function(Z, t){
        Xpi_Q(Z,t,lambda = lambda_Q, T_=T_, beta = beta, a=a, b=b, k=k, C=C, delta = delta_Q)
    }
    
    with_seed(seed, {
        ## sim dW / W /
        dW <- rnorm(Nseq, sd = sqrt(dt))
        W <- cumsum(dW)
        Z <- exp(-0.5*beta^2*tseq + beta*W)
        ## sim X_pi
        XU <- map2_dbl(Z, tseq, Xp_U)
        XC <- map2_dbl(Z, tseq, Xp_C)
        XP <- map2_dbl(Z, tseq, Xp_P)
        XE <- map2_dbl(Z, tseq, Xp_E)
        XQ <- map2_dbl(Z, tseq, Xp_Q)
        ##
        df <- tibble(
            t = c(0, tseq),
            Z = c(x,x+a*tseq+sigma*W) - (1-tilde)*c(0, (a-b)*tseq),
            XU = c(x, XU) - (1-tilde)*c(0, (a-b)*tseq),
            XC = c(x, XC) - (1-tilde)*c(0, (a-b)*tseq),
            XP = c(x, XP) - (1-tilde)*c(0, (a-b)*tseq),
            XE = c(x, XE) - (1-tilde)*c(0, (a-b)*tseq),
            XQ = c(x, XQ) - (1-tilde)*c(0, (a-b)*tseq),
        )
        df <- df |> pivot_longer(
            -1,
            names_to = "Constraints",
            values_to = "path"
        )
    })
    
    inputs <- c(a = a, b = b, sigma = sigma, T_ = T_, x = x, k0 = k0, C0 = C0, epsilon = epsilon, nu = nu)
    
    return (list(paths = df, seed = seed, tilde = tilde, inputs = inputs))
}

plot_paths_X <- function(paths_X) {
    # Input paths is an output of get_paths_X
    
    seed <- paths_X$seed
    path <- paths_X$paths
    tilde <- paths_X$tilde
    a <- paths_X$inputs[1]
    b <- paths_X$inputs[2]
    T_ <- paths_X$inputs[4]
    k0 <- paths_X$inputs[6]
    C0 <- paths_X$inputs[7]
    
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    gplot <- ggplot(path, aes(x = t, y = path, col = Constraints)) + 
        geom_line() + 
        geom_segment(aes(x=0,xend=T_, y=k, yend=tilde*k0 + (1-tilde)*k), linetype = "dashed", color = "black") +
        geom_segment(aes(x=0,xend=T_, y=C, yend=tilde*C0 + (1-tilde)*C), linetype = "dashed", color = "black") +
        annotate(geom = "text",
                 x = T_, y = tilde*k0 + (1-tilde)*k,
                 label = "italic(k)",parse=TRUE,
                 hjust = -0.5, vjust = 0.25) +
        annotate(geom = "text",
                 x = T_, y = tilde*C0 + (1-tilde)*C,
                 label = "italic(C)",parse=TRUE,
                 hjust = -0.5, vjust = 0.25) +
        labs(
            x = "Time",
            y = "Path",
            title = "Optimal portfolio value under different constraints",
            subtitle = sprintf("Seed: %s", seed)
        ) +
        theme_bw()
    
    return (gplot)
} 

# pi path ------------------------------------------------------------------------------------------------------------

calc_paths_pi <- function(seed, Nseq, 
                         a, b, sigma, T_, x, k0, C0, epsilon, nu,
                         tilde = T) {
    beta = -1*b/sigma
    ## 
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    ##
    
    tseq <- seq(0,T_,length.out = Nseq+1)[-1]
    tseq[length(tseq)] <- tseq[length(tseq)] - 0.00001
    dt <- tseq[2]-tseq[1]
    
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    pi_t_U <- function(Z, t){
        pi_U(Z,t,lambda = lambda_U, T_=T_, beta = beta, sigma=sigma)
    }
    
    lambda_C <- par_constrained(a, b, sigma, T_, x, k0, C0)$lambda_C
    pi_t_C <- function(Z, t){
        pi_C(Z,t,lambda = lambda_C, T_=T_, beta = beta, sigma=sigma, k=k, C=C)
    }
    
    par_P <- par_prob_constraint(a, b, sigma, T_, x, k0, C0, epsilon)
    lambda_P <- par_P$lambda_P
    g1 <- par_P$g1
    g2 <- par_P$g2
    pi_t_P <- function(Z, t){
        pi_P(Z,t,lambda = lambda_P, T_=T_, beta = beta, sigma=sigma, k=k, C=C, g1=g1, g2=g2)
    }
    
    par_E <- par_ES_constraint_P(a, b, sigma, T_, x, k0, C0, nu)
    lambda_E <- par_E$lambda_E
    h1 <- par_E$h1
    h2 <- par_E$h2
    pi_t_E <- function(Z, t){
        pi_E(Z,t,lambda = lambda_E, T_=T_, beta = beta, sigma=sigma, k=k, C=C, h1=h1, h2=h2)
    }
    
    par_Q <- par_ES_constraint_Q(a, b, sigma, T_, x, k0, C0, nu)
    lambda_Q <- par_Q$lambda_Q
    delta_Q <- par_Q$delta_Q
    pi_t_Q <- function(Z, t){
        pi_Q(Z,t,lambda = lambda_Q, T_=T_, beta = beta, sigma=sigma, k=k, C=C, delta = delta_Q)
    }
    
    with_seed(seed, {
        ## sim dW / W /
        dW <- rnorm(Nseq, sd = sqrt(dt))
        W <- cumsum(dW)
        Z <- exp(-0.5*beta^2*tseq + beta*W)
        ## sim X_pi
        piU <- map2_dbl(Z, tseq, pi_t_U)
        piC <- map2_dbl(Z, tseq, pi_t_C)
        piP <- map2_dbl(Z, tseq, pi_t_P)
        piE <- map2_dbl(Z, tseq, pi_t_E)
        piQ <- map2_dbl(Z, tseq, pi_t_Q)
        ##
        df <- tibble(
            t = c(0, tseq),
            pit_X = 0,
            pit_U = c(0, piU),
            pit_C = c(0, piC),
            pit_P = c(0, piP),
            pit_E = c(0, piE),
            pit_Q = c(0, piQ),
        )
        df <- df |> pivot_longer(
            -1,
            names_to = "Constraints",
            values_to = "path"
        )
    })
    
    inputs <- c(a = a, b = b, sigma = sigma, T_ = T_, x = x, k0 = k0, C0 = C0, epsilon = epsilon, nu = nu)
    
    return (list(paths = df, seed = seed, inputs = inputs))
}


plot_paths_pi <- function(paths_pi) {
    # Input paths is an output of get_paths_pi
    
    seed <- paths_pi$seed
    # pi will be 0 at time 0 for all strategies
    path <- paths_pi$paths%>% 
        dplyr::filter(t != 0)
    a <- paths_pi$inputs[1]
    b <- paths_pi$inputs[2]
    T_ <- paths_pi$inputs[4]
    k0 <- paths_pi$inputs[6]
    C0 <- paths_pi$inputs[7]
    
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    ylim_buffer = 0.001
    
    gplot <- ggplot(path, aes(x = t, y = path, col = Constraints)) + 
        geom_line() + 
        ylim(min(path$path) - ylim_buffer, 1 + ylim_buffer) +
        labs(
            x = "Time",
            y = "Path",
            title = "Optimal pi under different constraints",
            subtitle = sprintf("Seed: %s", seed)
        ) +
        geom_hline(yintercept = 1, lty = "dashed")
        theme_bw()
    
    return (gplot)
} 

# ---------------------------------------------------------------------------------------------------------------------------------
# output check for Pi ------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------

calc_paths_pi_check <- function(seed, Nseq, 
                                a, b, sigma, T_, x, k0, C0, epsilon, nu) {
    beta = -1*b/sigma
    ## 
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    ##
    
    tseq <- seq(0,T_,length.out = Nseq+1)[-1]
    tseq[length(tseq)] <- tseq[length(tseq)] - 0.00001
    dt <- tseq[2]-tseq[1]
    
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    pi_t_U <- function(Z, t){
        pi_U(Z,t,lambda = lambda_U, T_=T_, beta = beta, sigma=sigma)
    }
    
    lambda_C <- par_constrained(a, b, sigma, T_, x, k0, C0)$lambda_C
    pi_t_C <- function(Z, t){
        pi_C(Z,t,lambda = lambda_C, T_=T_, beta = beta, sigma=sigma, k=k, C=C)
    }
    
    par_P <- par_prob_constraint(a, b, sigma, T_, x, k0, C0, epsilon)
    lambda_P <- par_P$lambda_P
    g1 <- par_P$g1
    g2 <- par_P$g2
    pi_t_P <- function(Z, t){
        pi_P(Z,t,lambda = lambda_P, T_=T_, beta = beta, sigma=sigma, k=k, C=C, g1=g1, g2=g2)
    }
    
    par_E <- par_ES_constraint_P(a, b, sigma, T_, x, k0, C0, nu)
    lambda_E <- par_E$lambda_E
    h1 <- par_E$h1
    h2 <- par_E$h2
    pi_t_E <- function(Z, t){
        pi_E(Z,t,lambda = lambda_E, T_=T_, beta = beta, sigma=sigma, k=k, C=C, h1=h1, h2=h2)
    }
    
    par_Q <- par_ES_constraint_Q(a, b, sigma, T_, x, k0, C0, nu)
    lambda_Q <- par_Q$lambda_Q
    delta_Q <- par_Q$delta_Q
    pi_t_Q <- function(Z, t){
        pi_Q(Z,t,lambda = lambda_Q, T_=T_, beta = beta, sigma=sigma, k=k, C=C, delta = delta_Q)
    }
    
    with_seed(seed, {
        ## sim dW / W /
        dW <- rnorm(Nseq, sd = sqrt(dt))
        W <- cumsum(dW)
        Z <- exp(-0.5*beta^2*tseq + beta*W)
        ## sim X_pi
        piU <- map2_dbl(Z, tseq, pi_t_U)
        piC <- map2_dbl(Z, tseq, pi_t_C)
        piP <- map2_dbl(Z, tseq, pi_t_P)
        piE <- map2_dbl(Z, tseq, pi_t_E)
        piQ <- map2_dbl(Z, tseq, pi_t_Q)
        ##
        df <- tibble(
            t = c(0, tseq),
            Z = c(x, x + (a-b*0)*tseq + (1-0)*sigma*W),
            XU = c(x, x + (a-b*piU)*tseq + (1-piU)*sigma*W),
            XC = c(x, x + (a-b*piC)*tseq + (1-piC)*sigma*W),
            XP = c(x, x + (a-b*piP)*tseq + (1-piP)*sigma*W),
            XE = c(x, x + (a-b*piE)*tseq + (1-piE)*sigma*W),
            XQ = c(x, x + (a-b*piQ)*tseq + (1-piQ)*sigma*W),
        )
        df <- df |> pivot_longer(
            -1,
            names_to = "Constraints",
            values_to = "path"
        )
    })
    
    inputs <- c(a = a, b = b, sigma = sigma, T_ = T_, x = x, k0 = k0, C0 = C0, epsilon = epsilon, nu = nu)
    
    return (list(paths = df, seed = seed, inputs = inputs))
}

# ---------------------------------------------------------------------------------------------------------------------------------
# Plot both Xpi and pi
# ---------------------------------------------------------------------------------------------------------------------------------

plot_paths_pi_check <- function(paths_pi_check) {
    # Input paths is an output of get_paths_pi
    
    seed <- paths_pi_check$seed
    path <- paths_pi_check$paths
    a <- paths_pi_check$inputs[1]
    b <- paths_pi_check$inputs[2]
    T_ <- paths_pi_check$inputs[4]
    k0 <- paths_pi_check$inputs[6]
    C0 <- paths_pi_check$inputs[7]
    
    k <- k0 - (a-b)*T_
    C <- C0 - (a-b)*T_
    
    gplot <- ggplot(path, aes(x = t, y = path, col = Constraints)) + 
        geom_line() + 
        geom_segment(aes(x=0,xend=T_, y=k, yend=k0), linetype = "dashed", color = "black") +
        geom_segment(aes(x=0,xend=T_, y=C, yend=C0), linetype = "dashed", color = "black") +
        annotate(geom = "text",
                 x = T_, y = k0,
                 label = "italic(k)",parse=TRUE,
                 hjust = -0.5, vjust = 0.25) +
        annotate(geom = "text",
                 x = T_, y = C0,
                 label = "italic(C)",parse=TRUE,
                 hjust = -0.5, vjust = 0.25) +
        labs(
            x = "Time",
            y = "Path",
            title = "Optimal portfolio value under different constraints (Calculated using pi)",
            subtitle = sprintf("Seed: %s", seed)
        ) +
        theme_bw()
    
    return (gplot)
}

plot_paths_X_pi <- function(seed, Nseq, rel_heights = c(0.05, 1)){
    
    tseq <- seq(0,T_,length.out = Nseq+1)[-1]
    tseq[length(tseq)] <- tseq[length(tseq)] - 0.00001
    dt <- tseq[2]-tseq[1]
    
    lambda_U <- par_unconstrained(a, b, sigma, T_, x, k0, C0)$lambda_U
    Xp_U <- function(Z, t){
        Xpi_U(Z,t,lambda = lambda_U, T_=T_, beta = beta, a=a, b=b, k=k)
    }
    
    lambda_C <- par_constrained(a, b, sigma, T_, x, k0, C0)$lambda_C
    Xp_C <- function(Z, t){
        Xpi_C(Z,t,lambda = lambda_C, T_=T_, beta = beta, a=a, b=b, k=k, C=C)
    }
    
    par_P <- par_prob_constraint(a, b, sigma, T_, x, k0, C0, epsilon)
    lambda_P <- par_P$lambda_P
    c_P <- par_P$c_P
    g1 <- par_P$g1
    g2 <- par_P$g2
    Xp_P <- function(Z, t){
        Xpi_P(Z,t,lambda = lambda_P, T_=T_, beta = beta, a=a, b=b, k=k, C=C, g1=g1, g2=g2)
    }
    
    par_E <- par_ES_constraint_P(a, b, sigma, T_, x, k0, C0, nu)
    lambda_E <- par_E$lambda_E
    gamma_E <- par_E$gamma_E
    h1 <- par_E$h1
    h2 <- par_E$h2
    Xp_E <- function(Z, t){
        Xpi_E(Z,t,lambda = lambda_E, T_=T_, beta = beta, a=a, b=b, k=k, C=C, h1=h1, h2=h2)
    }
    
    par_Q <- par_ES_constraint_Q(a, b, sigma, T_, x, k0, C0, nu)
    lambda_Q <- par_Q$lambda_Q
    delta_Q <- par_Q$delta_Q
    Xp_Q <- function(Z, t){
        Xpi_Q(Z,t,lambda = lambda_Q, T_=T_, beta = beta, a=a, b=b, k=k, C=C, delta = delta_Q)
    }
    
    get_paths_df <- function(seed){
        with_seed(seed, {
            ## sim dW / W /
            dW <- rnorm(Nseq, sd = sqrt(dt))
            W <- cumsum(dW)
            Z <- exp(-0.5*beta^2*tseq + beta*W)
            ## sim X_pi
            XU <- map2_dbl(Z,tseq,Xp_U)
            XC <- map2_dbl(Z,tseq,Xp_C)
            XP <- map2_dbl(Z,tseq,Xp_P)
            XE <- map2_dbl(Z,tseq,Xp_E)
            XQ <- map2_dbl(Z,tseq,Xp_Q)
            ##
            Z <- exp(-0.5*beta^2*tseq + beta*W)
            df <- tibble(
                t = c(0, tseq),
                `Original X` = x + c(0,a*tseq+sigma*W),
                Z = c(1,Z),
                `Unconstrained` = c(x, XU),
                `Strict Constraint` = c(x, XC),
                `VaR Constraint` = c(x, XP),
                `P-ES Constraint` = c(x, XE),
                `Q-ES Constraint` = c(x, XQ),
                `pit_Original X` = 0,
                `pit_Unconstrained` = pi_U(Z,t,lambda_U,T_,beta,sigma),
                `pit_Strict Constraint` = pi_C(Z,t,lambda_C,T_,beta,sigma,k,C),
                `pit_VaR Constraint` = pi_P(Z,t,lambda_P,T_,beta,sigma,k,C,g1,g2),
                `pit_P-ES Constraint` = pi_E(Z, t,lambda_E, T_, beta, sigma, k, C, h1, h2),
                `pit_Q-ES Constraint` = pi_Q(Z, t,lambda_Q, T_, beta, sigma, k, C, delta_Q)
            )
            
        })
        df
    }
    
    # new version
    df <- get_paths_df(seed)
    # suff <- c("", "U", "C", "P", "E", "Q")
    # suff2 <- c("X", "U", "C", "P", "E", "Q")
    suff <- c("Original X", "Unconstrained", "Strict Constraint", "VaR Constraint", "P-ES Constraint", "Q-ES Constraint")
    suff2 <- c("Original X", "Unconstrained", "Strict Constraint", "VaR Constraint", "P-ES Constraint", "Q-ES Constraint")
    sp <- tibble(
        .name = c(
            paste0(suff),
            paste0("pit_",suff2)
        ),
        .value = c(suff2,suff2),
        type = c(
            rep("path", length(suff)),
            rep("proportion", length(suff))
        )
    )
    df_plot <- pivot_longer_spec(df,sp) %>%  
        pivot_longer(
            cols = all_of(suff2), names_to = "constraints", values_to = "y"
        )
    ##
    linesize = 0.6
    legend_linesize = 1.2
    # from scale::show_col(hue_pal()(6))
    colors = c("Original X" = "#F8766D", 
               "Unconstrained" = "#B79F00", 
               "Strict Constraint" = "#00BA38", 
               "VaR Constraint" = "#00BFC4", 
               "P-ES Constraint" = "#619CFF",
               "Q-ES Constraint" = "#F564E3")
    
    fig1 <- ggplot(df_plot %>% dplyr::filter(type == "path"), aes(x=t, y=y, color = constraints)) + 
        geom_line(size = linesize) + 
        geom_segment(aes(x=0,xend=T_, y=k,yend=k0), linetype = "dashed", color = "black") +
        geom_segment(aes(x=0,xend=T_, y=C,yend=C0), linetype = "dashed", color = "black") +
        annotate(geom = "text",
                 x = T_, y = k0,
                 label = "italic(k)",parse=TRUE,
                 hjust = -0.5, vjust = 0.25) +
        annotate(geom = "text",
                 x = T_, y = C0,
                 label = "italic(C)",parse=TRUE,
                 hjust = -0.5, vjust = 0.25) +
        labs(
            y = expression(paste("Surplus value ", italic(X)[t]))
        ) +
        theme_bw() + 
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank()) +
        guides(color = guide_legend(override.aes = list(size = legend_linesize))) +
        scale_colour_manual(name = "Constraints", values = colors)
    ##
    fig2 <- ggplot(df_plot %>% dplyr::filter(type == "proportion"), aes(x=t, y=y, color = constraints)) + 
        geom_line(size = linesize) + 
        geom_hline(yintercept = 1, linetype = "dashed") +
        labs(
            x = expression(paste("Time ", italic("t"))),
            y = expression(paste("Reinsurance control ", italic(pi)[t]))
        ) + 
        theme_bw() +
        guides(color = guide_legend(override.aes = list(size = legend_linesize))) +
        scale_colour_manual(name = "Constraints", values = colors)
    ## use cowplot
    plot_title <- ggdraw() + 
        draw_label(
           "Optimal portfolio value under different Constraints",
           fontface = 'bold', x = 0, hjust = 0
        )
    plot_subtitle <- ggdraw() + 
        draw_label(
            sprintf("Seed: %s", seed),
            size = 10,
            x = 0, hjust = 0
        )
    plot_main <- plot_grid(fig1,fig2, ncol = 1)
    ## return
    plot_grid(
        # plot_title,
        plot_subtitle,
        plot_main,
        ncol = 1,
        rel_heights = rel_heights
    )
}