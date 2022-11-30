sample_phylo <- function(gamma_sus, gamma_u, gamma_t, n_tip, traj_seed, phy_seed) {
    set.seed(traj_seed)

    usage_df <- simulated_usage()
    usage_f <- yearly_usg_stepfunc(usage_df$usage, usage_df$time)
    
    t0 <- min(usage_df$time)
    tmax <- max(usage_df$time)

    n_step <- (tmax-t0)/1e-4
    
    N_0 <- 3e6
    n_strains <- 3

    t0 <- min(usage_df$time)
    tmax <- max(usage_df$time)
    usage_f <- yearly_usg_stepfunc(usage_df$usage, usage_df$time)


    incidence <- 1.5e5 / N_0
    mu <- incidence+gamma_sus
    r0 <- mu/(gamma_sus)
    end_eq <- 1-1/r0
        
    N_f <- function(t) N_0*exp(log(2)*(t-t0)/(tmax-t0))

    t_change <- (t0 + floor(0.33*(tmax-t0)))

    beta_f <- function (t) if (t < t_change) 1.02*mu else 
        mu*(0.98 + 0.06*(t-t_change)/(tmax-t_change))

    I_0 <- rep(0, n_strains)
    while(!all(I_0 >= 100) || !(sum(I_0) < N_0)){
        Itmp <- floor(end_eq*N_0*rlnorm(n_strains-1, -4,.3))
        I_0 <- c(N_0*end_eq-sum(Itmp), Itmp)
    }

    S_0 <- N_0 - sum(I_0)
    i <- 1

    stopifnot(S_0 > 1)

    gamma_u_ext <- c(gamma_sus, gamma_sus, gamma_u)
    gamma_t_ext <- c(gamma_sus, gamma_sus, gamma_t)

    df <- simulate_epi_nstrain(t0, tmax, n_step, n_strains, S_0, I_0, usage_f, beta_f, N_f, gamma_u_ext, gamma_t_ext)
    while (!all(df[,paste0("I_", c(1:n_strains))] > 0)) {
        if (i > 10) stop("Stopping after 10 resampling attempts")
        i <- i+1
        df <- simulate_epi_nstrain(t0, tmax, n_step, n_strains, S_0, I_0, usage_f, beta_f, N_f, gamma_u_ext, gamma_t_ext)
    }

    set.seed(phy_seed)

    sampled_sus <- sample_phy_incidence(df, 2, n_tip, c(0:5), t0, tmax)
    sampled_res <- sample_phy_incidence(df, 3, n_tip, c(0:5), t0, tmax)

    phys <- list(sampled_sus$phy, sampled_res$phy)
    most_recent_samp <- list(sampled_sus$mr, sampled_res$mr)
 
    return(list(phys=phys, traj=df, most_recent_samp=most_recent_samp))
}

sample_phy_incidence <- function(df, strain_idx, n_tips, samp_yrs, t0, tmax) {
    I_col <- paste0("I_", strain_idx)
    inci_col <- paste0("births_", strain_idx)
    incidence_yr <- rev(sapply(c(t0:(tmax-1)), function(t) sum(df[which((df$t >= t) & (df$t < (t+1))), inci_col])))

    n_samp <- rmultinom(1,n_tips, incidence_yr[samp_yrs]/sum(incidence_yr[samp_yrs]))
    samp_times <- samp_yrs[which(n_samp>0)]
    n_samp <- n_samp[which(n_samp>0)]
    mr <- min(samp_times)

    phy <- simulate_phylos_nstrain(df, strain_idx, list(samp_times), list(n_samp))
    return(list(phy=phy[[1]],mr=mr))
}