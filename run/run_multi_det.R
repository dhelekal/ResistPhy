library(ggplot2)
library(reshape2)
library(viridis)
library(ResistPhy)
library(ape)
library(posterior)
library(cmdstanr)
library(bayesplot)
library(phylodyn)
library(latex2exp)

source("simu_phylo.R")

# test if there is at least one argument: if not, return an error
args <- commandArgs(TRUE)
print(args)
if (length(args)!=2) {
  stop("Need 2 arguments.", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = "out.txt"
}

idx_begin <- as.numeric(args[1])
extent <- as.numeric(args[2])

n_runs <- 50

run_idx <- c(idx_begin:(idx_begin + extent - 1))

usage_df <- simulated_usage()

set.seed(123456)
traj_seeds <- sample.int(n=100000, 50, F)
set.seed(457869)
phy_seeds <- sample.int(n=100000, 50, F)

q_u <- seq(from=1.0, to=1.2, length.out=n_runs)
q_t <- seq(from=1.0, to=0.5, length.out=n_runs)

t0 <- min(usage_df$time)
tmax <- max(usage_df$time)
gamma_sus <- 1/60.0*365.0
n_tip <- 200

gamma_u <- gamma_sus*q_u
gamma_t <- gamma_sus*q_t

phys <- list()
most_recent_samp <- list()

for (i in c(1:extent)) {
    idx <- run_idx[i]

    sim_out <- sample_phylo(gamma_sus, gamma_u[idx], gamma_t[idx], n_tip, traj_seeds[idx], phy_seeds[idx])
    sim_df <- sim_out$traj

    phys[[i]] <- sim_out$phys
    most_recent_samp[[i]] <- sim_out$most_recent_samp
}

n_warmup <- 2000
n_it <- 2000
recovery_data <- data.frame(variable=c(), value=c(), para_idx=c())
converged <- c()

for (i in c(1:extent)) {
    out <- infer_costs2(phys[[i]], 
                most_recent_samp[[i]],
                usage_df$usage,
                usage_df$time,
                t0,
                tmax,
                gamma_sus/365,
                365,
                n_iter=n_it, 
                n_warmup=n_warmup,
                model="deterministic",
                K=60,
                L=6.5,
                gamma_log_sd=0.1,
                seed=run_idx[i],
                stan_control=list(adapt_delta=.99,
                    max_treedepth=13,
                    parallel_chains=4,
                    chains=4,
                    refresh=0
                ))
    
    temp_df <- out$draws_df[,c("q_u[1]", "q_t[1]")]
    colnames(temp_df) <- c("q_u","q_t")

    temp_df <- melt(temp_df, measure.vars=c("q_u", "q_t"), 
                variable.name="variable", value.name="value")
    temp_df$para_idx <- run_idx[i]
    recovery_data <- rbind(recovery_data, temp_df)
    converged <- c(converged, out$converged)
}

saveRDS(recovery_data, paste0("det_para_",idx_begin,"_",extent, ".rds"))
saveRDS(converged, paste0("det_conv_",idx_begin,"_",extent, ".rds"))