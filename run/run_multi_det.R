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
if (length(args)!=3) {
  stop("Need 3 arguments.", call.=FALSE)
}

idx_begin <- as.numeric(args[1])
extent <- as.numeric(args[2])
mod_c <- as.numeric(args[3])

if(mod_c==1) {
  mod <- "nodecay"
} else if (mod_c==2) {
  mod <- "decay"
} else {
  stop("Invalid model choice", call.=FALSE)
}


n_runs <- 50
usage_df <- simulated_usage()

n_warmup <- 2000
n_it <- 2000
recovery_data <- data.frame(variable=c(), value=c(), para_idx=c())
converged <- c()

sim<-sim_n_traj(idx_begin, extent, n_runs, plot_traj=FALSE)

for (i in c(1:extent)) {
    out <- infer_costs2(sim$phys[[i]], 
                sim$most_recent_samp[[i]],
                usage_df$usage,
                usage_df$time,
                t0,
                tmax,
                gamma_sus/365,
                365,
                n_iter=n_it, 
                n_warmup=n_warmup,
                model=mod,
                K=60,
                L=6.5,
                gamma_log_sd=0.1,
                seed=sim$run_idx[i],
                stan_control=list(adapt_delta=.99,
                    max_treedepth=13,
                    parallel_chains=4,
                    chains=4,
                    refresh=0
                ))
      if (mod == "nodecay"){
        dnames <- c("q_u[1]", "q_t[1]")
        vnames <- c("q_u","q_t")
      } else {
        dnames <- c("q_u[1]", "q_t[1]", "phi[1]", "p[1]")
        vnames <- c("q_u","q_t", "phi", "p")
      }
      temp_df <- out$draws_df[,dnames]
      colnames(temp_df) <- vnames

      temp_df <- melt(temp_df, measure.vars=vnames, 
        variable.name="variable", value.name="value")

      temp_df$para_idx <- sim$run_idx[i]
      recovery_data <- rbind(recovery_data, temp_df)
      converged <- c(converged, out$converged)
}

saveRDS(recovery_data, paste0("det_para_",idx_begin,"_",extent, ".rds"))
saveRDS(converged, paste0("det_conv_",idx_begin,"_",extent, ".rds"))
