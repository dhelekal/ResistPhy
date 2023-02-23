costs_fit2 <- function(fit, 
    n_lineages, 
    time_vecs, 
    phys, 
    idx_begin, 
    idx_end, 
    time_scale, 
    gamma_guess, 
    gamma_log_sd, 
    n_iter, 
    n_warmup) {

    draws_a <- as_draws_array(fit$draws())
    ddf <- as.data.frame(as_draws_df(fit$draws()))
    np_cp <- nuts_params(fit)

    sampler_diagnostics <- fit$sampler_diagnostics()
    sampler_diagnostics <- as_draws_df(sampler_diagnostics)

    fit_summary <- fit$summary(c("alpha", 
                        "rho", 
                        "f_tilde", 
                        "gamma_sus_tilde",
                        "q_tilde", 
                        "I_0_hat"))

    converged <- TRUE
    ess_bulk <- fit_summary$ess_bulk

    if(!all(fit_summary$rhat < 1.05)) {
        warning("Poor convergence detected: unsatisfactory rhat")
        converged <- FALSE
    }
    if (!all(ess_bulk > 500)) {
        warning("Poor convergence detected: Bulk ESS lower than 500")
        converged <- FALSE
    }
    if(length(which(sampler_diagnostics$divergent__ > 1e-8))>1) {
        warning(sprintf("%d/%d divergent transitions detected",
                length(which(sampler_diagnostics$divergent__ > 1e-8)),
                length(sampler_diagnostics$divergent__)))
        converged <- FALSE
    }

    out <- list(
        draws_arr=draws_a,
        draws_df=ddf,
        n_lineages=n_lineages,
        time_vecs=time_vecs,
        phys=phys,
        idx_begin=idx_begin,
        idx_end=idx_end,
        gamma_log_sd=gamma_log_sd,
        gamma_guess=gamma_guess,
        time_scale=time_scale,
        n_iter=n_iter,
        n_warmup=n_warmup,
        stan_fit=fit,

        converged=converged
    )

    attr(out, "class") <- "costsFit2"
    return(out)
}

#' Print fit summary
#' @param o costsFit2 object
#' @export
print.costsFit2 <- function (o, ...) {
    cat(paste("\nResistPhy NUTS posterior fit\n\n"))  
    cat(paste("\nFields: ",  names(o), "\n"))
    cat(paste("\nNumber of strains: ", o$n_lineages))
    cat(paste("\nNumber of iterations: ", o$n_iter))
    cat(paste("\nNumber of warmup iterations: ", o$n_warmup,"\n"))
    cat(paste("\nConvergence: ", o$converged,"\n"))
}

#' Plot traces for a subset of parameters using bayesplot
#' @param o costsFit2 object
#' @export
plot_traces <- function(o, ...) {
    mcmc_trace(o$draws_arr, pars = c("alpha", "rho", "gamma_sus_tilde"), regex_pars=c("I_0_hat.*", "q_tilde.*"),
           facet_args = list(ncol = 1, strip.position = "left"))
}

#' Plot posterior marginals for GP hyperparameters
#' @param o costsFit2 object
#' @export
plot_hyperpar_pairs <- function(o, ...) {
    blank <- plot_spacer()
    alpha_hist <- ggplot(o$draws_df, aes(x=alpha)) + 
                        geom_histogram(aes(y = ..density..), position="identity", bins=50, fill="lightsteelblue3") + 
                        geom_function(fun = function(x) dgamma(x, shape=4, rate=4), linetype="longdash", color="steelblue4") + 
                        labs(x="", y="", title=TeX("$\\alpha$")) +
                        theme_minimal() + thm1()
    rho_hist <- ggplot(o$draws_df, aes(x=rho)) + 
                        geom_histogram(aes(y = ..density..), position="identity", bins=50, fill="lightsteelblue3") + 
                        geom_function(fun = function(x) dgamma(1/x, shape=4.63, rate=2.21)/(x^2), linetype="longdash", color="steelblue4") + 
                        labs(x="", y="", title=TeX("$\\rho$")) +
                        theme_minimal() + thm1()
    
    p1 <- ggplot(o$draws_df, aes(x=alpha, y=rho)) + 
            geom_hex(bins = 20) + 
            scale_fill_viridis() + 
            theme_minimal() +
            thm2()

    alpha_hist + blank +
        p1 + rho_hist +
        plot_layout(ncol = 2)
}

#' Plot posterior marginals and scatter plots for recovery rates using bayesplot
#' @param res_lineage_idx index of resistant lineage to display
#' @param o costsFit2 object
#' @export
plot_epi_pairs <- function(o, res_lineage_idx, ...) {

    blank <- plot_spacer()

    q_u_name <- paste0("q_u[",res_lineage_idx,"]")
    q_t_name <- paste0("q_t[",res_lineage_idx,"]")

    gamma_sus_hist <- ggplot(o$draws_df, aes_(x=as.name("gamma_sus"))) + 
        geom_histogram(aes(y = ..density..), position="identity", bins=50, fill="lightsteelblue3") + 
        geom_function(fun = function(x) dlnorm(x, meanlog = log(o$gamma_guess), sdlog = o$gamma_log_sd), linetype="longdash", color="steelblue4") + 
        labs(x="", y="", title=TeX("$\\gamma$")) +
        theme_minimal() + thm1()


    q_u_hist <- ggplot(o$draws_df, aes_(x=as.name(q_u_name))) + 
                        geom_histogram(aes(y = ..density..), position="identity", bins=50, fill="lightsteelblue3") + 
                        geom_function(fun = function(x) dlnorm(x, meanlog = 0, sdlog = 0.5), linetype="longdash", color="steelblue4") + 
                        labs(x="", y="", title=TeX(paste0("$",q_u_name,"$"))) +
                        theme_minimal() + thm1()
    q_t_hist <- ggplot(o$draws_df, aes_(x=as.name(q_t_name))) + 
                        geom_histogram(aes(y = ..density..), position="identity", bins=50, fill="lightsteelblue3") + 
                        geom_function(fun = function(x) dlnorm(x, meanlog = 0, sdlog = 0.5), linetype="longdash", color="steelblue4") + 
                        labs(x="", y="", title=TeX(paste0("$",q_t_name,"$"))) +
                        theme_minimal() + thm1()

    p1 <- ggplot(o$draws_df, aes_(x=as.name("gamma_sus"), y=as.name(q_u_name))) + 
            geom_hex(bins = 20) + 
            scale_fill_viridis() + 
            theme_minimal() +
            thm2()
    p2 <- ggplot(o$draws_df, aes_(x=as.name("gamma_sus"), y=as.name(q_t_name))) + 
            geom_hex(bins = 20) + 
            scale_fill_viridis() + 
            theme_minimal() +
            thm2()
    p3 <- ggplot(o$draws_df, aes_(x=as.name(q_u_name), y=as.name(q_t_name))) + 
            geom_hex(bins = 20) + 
            scale_fill_viridis() + 
            theme_minimal() +
            thm2() 
    
    gamma_sus_hist + blank + blank +
    p1 + q_u_hist + blank +
    p2 + p3 + q_t_hist +
    plot_layout(ncol = 3)
}

#' Plot q parameter marginals for all resistant strains
#' @param o costsFit2 object
#' @export
plot_qpairs <- function(o, ...) {
    blank <- plot_spacer()

    ddf <- o$draws_df
    n_lineages <- o$n_lineages

    lin_data <- list()

    for (i in c(1:(n_lineages-1))) {
        lin_names <- paste0("q_", c("u","t"),"[",i,"]")
        lin_df <- ddf[,lin_names]
        colnames(lin_df) <- c("q_u","q_t")
        lin_df$lineage <- i
        lin_data[[i]] <- lin_df
    }

    plot_df <- do.call(rbind, lin_data)
    plot_df$lineage <- factor(plot_df$lineage) 
    h1 <- ggplot(plot_df, aes(x=q_u, color=lineage, fill=lineage, group=lineage)) + 
        geom_density(alpha=0.4,lwd=0.8) + 
        labs(x="", y="", title=TeX(paste0("$","q_U","$"))) +
        scale_color_brewer(palette="Paired") +
        scale_fill_brewer(palette="Paired") +
        theme_minimal() + (thm1() %+replace% theme(legend.position="none"))
    h2 <- ggplot(plot_df, aes(x=q_t, color=lineage, fill=lineage, group=lineage)) + 
        geom_density(alpha=0.4, lwd=0.8) + 
        scale_color_brewer(palette="Paired") +
        scale_fill_brewer(palette="Paired") +
        labs(x="", y="", fill = "Resistant Lineage", color="Resistant Lineage", title=TeX(paste0("$","q_T","$"))) +
        theme_minimal() + thm1()

    p <- ggplot(plot_df, aes(x=q_u, q_t, color=lineage, fill=lineage, group=lineage)) +
        geom_density_2d(bins=10, alpha=.5) +
        scale_color_brewer(palette="Paired") + 
        theme_minimal() + thm2()

    h1 + blank +
    p + h2 +
    plot_layout(ncol = 2)
}

#' Plot epidemic trajectory credible intervals and posterior draws
#' @param o costsFit2 object
#' @param lineage_index index of lineage to display, '1' indicates susceptible
#' @param n_draws Number of posterior epidemic trajectory draws to sample
#' @export
plot_epi_dynamic <- function(o, lineage_index, n_draws=40, ...) {
    ts <- o$time_vecs[[lineage_index]]
    
    lineage_indices <- c(o$idx_begin[lineage_index]:o$idx_end[lineage_index])
    I_names <- sapply(lineage_indices,function(i) paste0("traj_vec[",i,"]"))
    ytitle <- "Prevalence"
    plot_timeseries(o, lineage_index, ts, I_names, ytitle, n_draws)
}

#' Plot Ne trajectories credible intervals and posterior draws
#' @param o costsFit2 object
#' @param lineage_index index of lineage to display, '1' indicates susceptible
#' @param n_draws Number of posterior Ne trajectory draws to sample
#' @export
plot_Ne <- function(o,  lineage_index, n_draws=40, ...) {
    ts <- o$time_vecs[[lineage_index]]
    
    lineage_indices <- c(o$idx_begin[lineage_index]:o$idx_end[lineage_index])
    I_names <- sapply(lineage_indices,function(i) paste0("Ne[",i,"]"))
    ytitle <- "Ne(t)"
    plot_timeseries(o, lineage_index, ts, I_names, ytitle, n_draws)
}

#' Plot growth rate r(t) credible intervals and posterior draws
#' @param o costsFit2 object
#' @param lineage_index index of lineage to display, '1' indicates susceptible
#' @param n_draws Number of posterior r(t) draws to sample
#' @export
plot_rt <- function(o, lineage_index, n_draws=40, ...) {
    ts <- o$time_vecs[[lineage_index]]
    
    lineage_indices <- c(o$idx_begin[lineage_index]:o$idx_end[lineage_index])
    I_names <- sapply(lineage_indices,function(i) paste0("r_t[",i,"]"))
    ytitle <- "r(t)"
    plot_timeseries(o, lineage_index, ts, I_names, ytitle, n_draws)
}

#' Plot Rt_res / Rt_sus ratio posterior probability curve for a given threshold
#' @description f(U) = P[rt_sus - rt_res > c | u(t) = U]
#' @param o costsFit2 object
#' @param res_lineage_idx index of resistant lineage to display
#' @param rr_threshold r_t ratio for which to plot the curve
#' @param n_breaks number of usage discretisation breaks
#' @export
plot_rr_curve <- function(o, res_lineage_idx, rr_threshold=0.0, n_breaks = 100) {
    q_names <- paste0(c("q_u[", "q_t["), res_lineage_idx, "]")
    q_df <- o$draws_df[,q_names]
    n_samp <- nrow(q_df)
    usage_thresholds <-  apply(q_df, 1, function(x) (rr_threshold - x[1])/(x[2]-x[1]))
    usage_breaks <- seq(from=0, to=1, length.out=n_breaks)
    probs <- sapply(usage_breaks, function(x) length(which(usage_thresholds > x))/n_samp)

    p_df <- data.frame(usage=usage_breaks, prob=probs)

    ggplot(p_df) + 
    geom_step(aes(x=usage, y=prob), fill="lightsteelblue3", color="lightsteelblue3") +
    labs(x="Usage", y="Posterior Probability", title=TeX(sprintf("Posterior Probability of $r_r(t) - r_s(t) < %.2f$", rr_threshold))) +
    theme_minimal() +
    thm3(aspect.ratio=1) 
}

#' Plot rt_sus - rt_res difference posterior probability heatmap
#' @description f(U, c) = P[rt_sus - rt_res > c | u(t) = U]
#' @param o costsFit2 object
#' @param res_lineage_idx index of resistant lineage to display
#' @param n_breaks number of usage discretisation breaks
#' @param min_disp_prob minimum probability display cutoff
#' @export
plot_rr_map <- function(o, res_lineage_idx, n_breaks = 100, min_disp_prob=0.5) {
    q_names <- paste0(c("q_u[", "q_t["), res_lineage_idx, "]")
    q_df <- o$draws_df[,q_names]
    n_samp <- nrow(q_df)
    usage_breaks <- seq(from=0, to=1, length.out=n_breaks)
    max_cost <- max(q_df$q_u)

    if (max_cost <= 0) {
         warning("Posterior difference between resistant and susceptible growth rates r(t) never drops below 0.
            This suggests no transmission cost to resistance")
    }

    rr_breaks <- seq(from=-max_cost, to = 0, length.out = n_breaks)
    p_df <- data.frame(prob=c(), usage=c(), rr=c())

    for (rr in rr_breaks) {
        usage_thresholds <-  apply(q_df, 1, function(x) (c - x[1])/(x[2]-x[1]))
        probs <- sapply(usage_breaks, function(x) length(which(usage_thresholds > x))/n_samp)
        tmp_df <- data.frame(prob = probs, usage=usage_breaks, rr=rr)
        p_df <- rbind(p_df, tmp_df)
    }

    x_upper <- max(p_df$usage[which(p_df$prob >= min_disp_prob)]) 
    y_lower <- min(p_df$rr[which(p_df$prob >= min_disp_prob)])

    breaks <- c(0, 0.50, 0.75, 0.90, 0.95, 0.99, 0.999, 1)

    ggplot(p_df) + 
    geom_contour_filled(aes(x=usage, y=rr, z=prob), breaks=breaks) +
    #geom_tile(aes(x=usage, y=rr, fill=prob)) +
    scale_color_viridis() +
    labs(x="Usage", y="C", fill="Posterior Probability", title=TeX(sprintf("Posterior Probability of $r_r(t) - r_s(t) < c$ given usage level"))) +
    xlim(0, x_upper) +
    ylim(y_lower, 1) +
    theme_minimal() +
    thm3(aspect.ratio=1) 
}

#' Plot posterior predictive draws and data overlay for the lineages through time function
#' @param o costsFit2 object
#' @param lineage_index index of lineage to display, '1' indicates susceptible
#' @param n_draws number of trajectories to simulate and overlay
#' @export
plot_ppcheck_At <- function(o, lineage_index, n_draws=100, ...) {
    ts <- o$time_vecs[[lineage_index]]
    t_max <- max(ts)-min(ts)
    n <- length(ts)

    phy <- o$phys[[lineage_index]]

    lineage_indices <- c(o$idx_begin[lineage_index]:o$idx_end[lineage_index])
    Ne_names <- sapply(lineage_indices,function(i) paste0("Ne[",i,"]"))

    row_count <- nrow(o$draws_df)
    stopifnot("n_draws must not exceed the number of MCMC draws"=n_draws <= row_count)
    draws <- sample.int(row_count, n_draws)

    pp_df <- draw_At_samples(phy, draws, o$draws_df, Ne_names, ts, t_max, 1.0)

    data_df <- compute_At(phy$samp_times, phy$n_samp, phy$coal_times, t_max)

    ggplot(pp_df) + geom_step(aes(x=t, y=At, group=draw), alpha=0.1, color="gray", size=0.2) + 
                    geom_step(data=data_df, aes(x=t, y=At), color="black", size=0.7, linetype="dashed") +
                    labs(x="Time", y="A(t)") + 
                    scale_y_log10() + 
                    theme_minimal() + 
                    thm3()

}

draw_At_samples <- function(phy, draws, draws_df, Ne_names, eval_ts, t_max, time_scale) {
    pp_df <- data.frame(t=c(), At=c(), draw=c())
    n_draws <- length(draws)
    for (i in c(1:n_draws)) {
        draw_idx <- draws[i]
        sample_phy <- sim_coal(phy$samp_times, phy$n_samp, unlist(draws_df[draw_idx,Ne_names])/time_scale, eval_ts)
        temp_df <- compute_At(sample_phy$samp_times, sample_phy$n_samp, sample_phy$coal_times,
                    t_max)
        temp_df$draw <- i;
        pp_df <- rbind(pp_df,temp_df)
        
    }
    return(pp_df)
}


plot_timeseries <- function(o, lineage_index, times, varnames, ytitle, n_draws) {
    ts <- times
    n <- length(ts)
    
    ddf <- o$draws_df
    
    # Compute CIs
    y_lo <- sapply(varnames, function (i) quantile(ddf[[i]],c(0.025)))
    y_hi <- sapply(varnames, function (i) quantile(ddf[[i]],c(0.975)))
    y_m <- sapply(varnames, function (i) median(ddf[[i]]))

    p_df <- data.frame(t=ts, 
                y=y_m,
                y_hi=y_lo, 
                y_lo=y_hi)
    
    # Compute posterior draws
    draw_idx <- sample.int(nrow(ddf), n_draws)
    names(ts) <- unlist(regmatches(varnames, regexec("\\d+",varnames)))
    
    traj_draws <- ddf[draw_idx, varnames]
    traj_draws$draw <- if(n_draws > 0) c(1:n_draws) else c()
    traj_draws <- melt(traj_draws, "draw", variable.name = "ID", value.name="y")
    draw_idx <- as.character(unlist(regmatches(traj_draws$ID, regexec("\\d+",traj_draws$ID))))
    traj_draws$t <- sapply(draw_idx, function(i) ts[[i]])

    ggplot(p_df)+geom_ribbon(aes(x=t, ymin=y_lo, ymax=y_hi), fill="lightgray", alpha=0.4) +
        geom_line(data=traj_draws, aes(x=t, y=y, color=draw, group=factor(draw)), alpha=0.4, size=0.1) +
        geom_line(aes(x=t, y=y), color="black", linetype="longdash") +
        scale_color_viridis() + 
        labs(title = paste("Lineage ", lineage_index), y=ytitle, x="Time") + 
        theme_minimal() + (thm3() %+replace% theme(legend.position="none"))
}

