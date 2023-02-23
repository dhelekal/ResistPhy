library(ggplot2)
library(reshape2)
library(viridis)
library(ResistPhy)
library(ape)
library(latex2exp)

source("simu_phylo.R")

#Simulate phylogenies and trajectories for a given subset of nruns simulations
sim_n_traj <- function (idx_begin, extent, n_runs, plot_traj=FALSE){
    out <- list()

    run_idx <- c(idx_begin:(idx_begin + extent - 1))

    usage_df <- simulated_usage()
    set.seed(123456)
    traj_seeds <- sample.int(n=100000, n_runs, F)
    set.seed(457869)
    phy_seeds <- sample.int(n=100000, n_runs, F)

    q_u <- seq(from=0.0, to = 1.3, length.out=n_runs)
    q_t <- seq(from=0.0, to = -2.9, length.out=n_runs)

    t0 <- min(usage_df$time)
    tmax <- max(usage_df$time)
    gamma_sus <- 1/60.0*365.0
    n_tip <- 200

    gamma_u <- gamma_sus+q_u
    gamma_t <- gamma_sus+q_t

    phys <- list()
    most_recent_samp <- list()

    if (plot_traj) {
        traj_df <- data.frame(t=c(), y=c(), strain=c(), para_idx=c())
        facet_labs <- sapply(c(1:n_runs), function(i) TeX(sprintf("$q_{u}$: %.3f $q_{t}$: %.3f", q_u[i], q_t[i])))
    }

    for (i in c(1:extent)) {
        idx <- run_idx[i]

        sim_out <- sample_phylo(gamma_sus, gamma_u[idx], gamma_t[idx], n_tip, traj_seeds[idx], phy_seeds[idx])

        phys[[i]] <- sim_out$phys
        most_recent_samp[[i]] <- sim_out$most_recent_samp

        if (plot_traj) #Downsample and retain trajectories
        {        
            sim_df <- sim_out$traj
            n <- nrow(sim_df)
            sub_idx <- seq.int(from=1, to=n, length.out=1e3)
            sim_df <- sim_df[sub_idx,]
            sim_df$Ne_s <- sim_df$I_2/2/sim_df$inf_rate 
            sim_df$Ne_r <- sim_df$I_3/2/sim_df$inf_rate 
            
            append_df <- melt(sim_df[,c("t","I_2","I_3")], "t", variable.name="strain", value.name="y") 

            append_df$para_idx <- i
            traj_df <- rbind(traj_df, append_df)
        }
    }

    out$run_idx <- run_idx
    out$phys <- phys
    out$most_recent_samp <- most_recent_samp
    out$q_u <- q_u
    out$q_t <- q_t

    if(plot_traj) {
        tplt <- ggplot(traj_df) + 
            geom_line(aes(x=t, y=y, color=strain)) + 
            scale_color_manual(values = c(I_3 = "red", I_2 = "blue"), labels=c("Resistant", "Susceptible")) + 
            facet_wrap(~para_idx, nrow=10, labeller=label_parsed) +
            labs(x="Time", y="Population Size") +
            theme_bw()
        out$tplt <- tplt
    }
    
    return(out)
}