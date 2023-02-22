#' Run fitness cost inference
#' @param phys A list of phylogenies. First phylogeny is assumed to susceptible
#' @param mrsts A list of most recent sampling times (in coalescent time) relative to t_begin 
#' @param ABX_usage AMR usage values as a proportion of primary treatment
#' @param ABX_times AMR usage time points
#' @param t_begin Analysis start date in natural time
#' @param t_end Analysis end date in natural time
#' @param gamma_guess Average recovery rate per day
#' @param time_scale_d The scale of one time unit for phylogenies and AMR data, in days
#' @param gamma_log_sd Log sd for strongly informative log-normal prior on average recovery rates
#' @param n_iter Number of sampling iterations
#' @param n_warmup Number of warmup iterations#
#' @param L HSGP approximation boundary factor
#' @param K HSGP approximation basis function count
#' @param seed RNG seed used by sampler
#' @param max_dt Maximum interval width for likelihood approximation on [-1,1] scale
#' @param stan_control List of additional arguments to pass to cmdstan. Accepted keywords: 'adapt_delta', 'max_treedepth', 'refresh', 'init', 'chains', 'parallel_chains'. See stan_defaults() for default values.
#' @export
#'
infer_costs2 <- function(phys, 
                        mrsts,
                        ABX_usage,
                        ABX_times,
                        t_begin,
                        t_end,
                        gamma_guess,
                        time_scale_d, 
                        gamma_log_sd=0.1,
                        n_iter = 2000, 
                        n_warmup = 2000, 
                        model="deterministic",
                        L = 6.5,
                        K = 65,
                        seed = 1,
                        max_dt=1e-3,
                        stan_control=stan_defaults()) {
    old_int <- c(t_begin, t_end)
    new_int <- c(-1,1)
    ABX_times <- rescale_times(ABX_times, old_int, new_int)
    scale_factor <- (old_int[2]-old_int[1])/2.0
    time_scale_new <- time_scale_d*scale_factor

    n_lineages <- length(phys)
    stopifnot(length(mrsts)==n_lineages)
    
    obs_total <- 0
    coal_total <- 0
    
    time_vecs <- list()
    obs_counts <- c()
    coal_counts <- c()

    all_times <- c()
    all_coal_idx <- c()
    all_w_t <- c()
    all_lineage_counts <- c()
    
    idx_begin <- c()
    idx_end <- c()

    obs_pos <- 1

    for (i in c(1:n_lineages)) {
        poi <- phy2poi(phys[[i]], mrsts[[i]], t_begin, t_end, max_dt)
        
        coal_idx_tmp <- which(poi$poi_counts > .5)
        times_tmp <- poi$process_times
        
        n_coal <- length(coal_idx_tmp)
        n_obs <- length(times_tmp)
        
        w_t_tmp <- poi$holding_times
        lineage_counts_tmp <- poi$lineage_counts

        all_times <- c(all_times, times_tmp)
        all_w_t <- c(all_w_t, w_t_tmp)
        all_lineage_counts <- c(all_lineage_counts, lineage_counts_tmp) 
        all_coal_idx <- c(all_coal_idx,coal_idx_tmp)

        idx_begin <- c(idx_begin, obs_pos)
        idx_end <- c(idx_end, obs_pos+n_obs-1)

        coal_counts <- c(coal_counts, n_coal)
        obs_counts <- c(obs_counts, n_obs)

        obs_total <- obs_total + n_obs
        coal_total <- coal_total + n_coal 

        obs_pos <- obs_pos + n_obs
        time_vecs[[i]] <- rescale_times(times_tmp, new_int, old_int)
    }

    data_list <- list(
        M=length(ABX_times),
        N_lineages=n_lineages,
        obs_total=obs_total,
        coal_total=coal_total,
        obs_counts=obs_counts,
        coal_counts=coal_counts,
        times=all_times,
        w_t=all_w_t,
        coal_index=all_coal_idx,
        lineage_counts=all_lineage_counts,
        usage_ts=ABX_times, 
        ABX_usage=ABX_usage, 
        time_scale=time_scale_new,
        unit_scale=time_scale_d,
        gamma_guess=gamma_guess,
        gamma_log_sd=gamma_log_sd,
        K=K,
        L=L
    )

    include_dir <- system.file('stan',package='ResistPhy',mustWork=T)

    if(model=="deterministic") {
        f <- system.file('stan',    
                        'n_strain_model_v8_deterministic.stan',
                        package='ResistPhy',
                        mustWork = T)
    } else {
        stop("Invalid Model Choice")
    }
    mod <- cmdstan_model(f, include_paths=include_dir, stanc_options = list("O1")) 

    fit <- mod$sample(data=data_list,
                    adapt_delta=stan_control$adapt_delta,
                    max_treedepth=stan_control$max_treedepth,
                    chains=stan_control$chains,
                    parallel_chains=stan_control$parallel_chains,
                    refresh=stan_control$refresh,
                    init=stan_control$init,
                    iter_warmup=n_warmup,
                    iter_sampling=n_iter,
                    seed=seed)

    phys_extracted <- list()

    for (i in c(1:n_lineages)) {
        phys_extracted[[i]] <- offset_extracted(extract_phy(phys[[i]]), mrsts[[i]])
    }

    out <- costs_fit2(fit, 
                    n_lineages,
                    time_vecs, 
                    phys_extracted, 
                    idx_begin,
                    idx_end,
                    time_scale_d, 
                    gamma_guess, 
                    gamma_log_sd, 
                    n_iter, 
                    n_warmup)
    return(out)
}


phy2poi <- function(phy, mrst, t_begin, t_end, max_dt=1e-3) {
    # Convert ape phylo to lists of:
    # - coalescent times
    # - sampling times 
    # - sampling counts
    extracted <- extract_phy(phy) 
    height_last_sam <- max(extracted$coal_times) - max(extracted$samp_times)
    # Offset times by most recent sampling time
    offset_phy <- offset_extracted(extracted, mrst)
    # Convert phylogeny to a list of event times, a list of event types, and list of lineage counts
    processed <- process_events(offset_phy$samp_times, offset_phy$n_samp, offset_phy$coal_times)
    # Truncate phylogeny to the usage data span
    # if MRCA exceeds the span
    span <- t_end-t_begin
    truncated <- truncate_phy(processed$event_times, processed$lineage_counts, processed$event_types, span)
    stopifnot(truncated$len > 2)
    # add P[No coalescence] (=sampling with 0 samples taken) event at t_begin (span in coalescent time)
    # with lineage count equal to lineage count before truncation point
    truncated$event_times <- c(truncated$event_times, span)
    truncated$lineages <- c(truncated$lineages, truncated$lineages_last)
    truncated$event_types <- c(truncated$event_types, 0)
    # reverse time to transform from coalescent time to natural time
    truncated$event_times <- t_end - truncated$event_times
    # rescale from t_begin to t_end interval to -1 to 1 interval
    old_int <- c(t_begin, t_end)
    new_int <- c(-1, 1)
    event_times_scaled <- rescale_times(truncated$event_times, old_int, new_int)
    stopifnot(max(abs(event_times_scaled[order(-event_times_scaled)] - event_times_scaled)) < 1e-8)
    # Reorder to natural time
    ord <- order(event_times_scaled)
    lineages <- truncated$lineages[ord]
    event_types <- truncated$event_types[ord]
    event_times <- event_times_scaled[ord]
    # Finally approximate phylogeny as a compound PPP, adding P[no event] intervals whenever inter-event times exceed max_dt
    n <- length(event_times_scaled)
    stopifnot(length(truncated$lineages)==n)
    stopifnot(length(truncated$event_types)==n)
    new_times <- c(event_times[1])
    new_lineages <- c(lineages[1])
    new_events <- c(event_types[1])
    stopifnot(event_types[1]==0)
    stopifnot(event_types[n]==0)
    stopifnot(lineages[n]==0)

    tol <- 1e-8

    for (i in c(2:n)) {
        next_t <- event_times[i]
        prev_t <- event_times[i-1]
        dt <- next_t - prev_t
        stopifnot(dt >= 0)
        if (dt > max_dt) {
            
            q <- (next_t-prev_t)/max_dt

            if (abs(trunc(q)-q) < tol) {
                n_step <- trunc(q)
            } else {
                n_step <- ceiling(q)
            }
            time_steps <- seq(from=prev_t, to=next_t, length.out=n_step+1)[-1]
            lins <- c(rep(lineages[i-1], n_step-1), lineages[i])
            evts <- c(rep(0, n_step-1),event_types[i])
            new_times <- c(new_times, time_steps)
            new_lineages <- c(new_lineages, lins)
            new_events <- c(new_events, evts)
        } else { 
            new_times <- c(new_times, event_times[i])
            new_lineages <- c(new_lineages, lineages[i])
            new_events <- c(new_events, event_types[i])
        }
    }

    m <- length(new_times)
    holding_times <- new_times[2:m] - new_times[1:(m-1)]
    process_times <- new_times[1:(m-1)]
    poi_counts <- new_events[1:(m-1)]
    lineage_counts <- new_lineages[1:(m-1)]

    return(list(holding_times=holding_times,
                process_times=process_times,
                poi_counts=poi_counts,
                lineage_counts=lineage_counts, 
                height_last_sam=height_last_sam))
}

rescale_times <- function(times, old_interval, new_interval){
    span_old <- old_interval[2]-old_interval[1]
    span_new <- new_interval[2]-new_interval[1]
    times_01 <- (times-old_interval[1])/span_old
    times_resc <- times_01*span_new + new_interval[1]
    stopifnot(all(times_resc <= new_interval[2]))
    stopifnot(all(times_resc >= new_interval[1]))
    return(times_resc)
}

#' Default stan arguments
#' @export
stan_defaults <- function() {
    return(list(adapt_delta=.99,
                max_treedepth=13,
                parallel_chains=4,
                chains=4,
                init=NULL))
}

offset_extracted <- function(extracted_phy, most_recent_sample) {
    # most recent sample as measured in time before present
    offset <- most_recent_sample
    return(list(coal_times=extracted_phy$coal_times+offset,
                samp_times=extracted_phy$samp_times+offset,
                n_samp=extracted_phy$n_samp))
}

extract_phy <- function(phy, merge_sam_threshold=1e-6) {
    phy_l <- makeNodeLabel(phy)
    labs <- c(phy_l$node.label, phy_l$tip.label)
    nodes <- nodeid(phy_l, labs)
    is_tip <- c(rep(FALSE, length(phy_l$node.label)), rep(TRUE, length(phy_l$tip.label)))

    times <- node.depth.edgelength(phy_l)
    times <- -times
    times <- times - min(times)
    times <- times[nodes]

    nodes.df <- data.frame(id=nodes, times=times, is_tip=is_tip, lab=labs)
    nodes.df <- nodes.df[order(nodes.df$times), ]

    coal_times <- c()
    samp_times <- c()
    n_samp <- c()
    i <- 1
    n <- nrow(nodes.df)

    while(i <= n) {
        time <- nodes.df$times[i]
        if (nodes.df$is_tip[i]) {
            k <- 1
            while(i < n && 
                nodes.df$is_tip[i+1] && 
                abs(time-nodes.df$times[i+1]) < merge_sam_threshold) {
                
                k <- k + 1
                i <- i + 1
            }
            i <- i + 1
            samp_times <- c(samp_times, time)
            n_samp <- c(n_samp, k)    
        } else {
            i <- i + 1
            coal_times <- c(coal_times, time)
        }
    }

    return(list(coal_times=coal_times,
                samp_times=samp_times,
                n_samp=n_samp))
}

truncate_phy <- function(event_times, lineages, event_types, span) {
    in_bounds <- which(event_times <= span)
    events_past_threshold <- which(event_times > span)
    if(length(events_past_threshold)>0) {
        lineages_last <- lineages[min(events_past_threshold)]
    } else {
        lineages_last <- 1
    }
    lineages_trunc <- lineages[in_bounds]
    event_types_trunc <- event_types[in_bounds]
    len <- length(lineages_trunc)
    times_trunc <- event_times[in_bounds]

    return(list(event_times=times_trunc, 
            lineages=lineages_trunc,
            event_types=event_types_trunc,
            len=len,
            lineages_last=lineages_last))
}

process_events <- function(samp_times, n_samp, coal_times) {
    stopifnot(samp_times == samp_times[order(samp_times)])
    stopifnot(coal_times == coal_times[order(coal_times)])

    event_types <- c(rep(0, length(samp_times)), rep(1, length(coal_times)))
    all_events <- c(samp_times, coal_times)
    lineage_counts <- rep(0, length(all_events))

    event_types <- event_types[order(all_events)]
    all_events <- all_events[order(all_events)]

    lin_count <- 0
    samp_idx <- 1

    for (i in c(1:length(all_events))) {
        lineage_counts[i] <- lin_count
        if (event_types[i] > 0.5) {
            ## coalescent event
            lin_count <- lin_count - 1
        } else {
            ## sampling event
            lin_count <- lin_count + n_samp[samp_idx]
            samp_idx <- samp_idx + 1
        }
    }
    return(list(event_times=all_events,
                event_types=event_types,
                lineage_counts=lineage_counts))
}