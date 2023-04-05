context("Coalescent")
test_that("Native poisson thinning method agrees with phylodyn",{
    set.seed(2)
    samp_times <- seq(from=0.0, to=2, by=0.1)
    n_samp <- 1+rbinom(length(samp_times), 20, 0.2)
    t0 <- 1
    tmax <- 20
    Ne_func <- function(t) exp(0.3*t) * 0.1
    Ne_times <- seq(from=t0, to=tmax,length.out=1000)
    Ne_vals <- sapply(Ne_times, Ne_func)

    dt <- Ne_times[2]-Ne_times[1]
    ne_step_func <- rev_time_step_func(Ne_vals, t0, tmax, dt)

    set.seed(2)
    res_phyd <- coalsim(samp_times, n_samp, ne_step_func, method="thin", lower_bound = min(Ne_vals))
    set.seed(2)
    res_native <- sim_coal(samp_times, n_samp, Ne_vals, Ne_times) 

    expect_equal(res_phyd$coal_times, res_native$coal_times)
})

context("Utils")
test_that("Phylogeny exttraction and offsetting agrees on simulated data",
    {   
        set.seed(1)
        samp_times <- seq(from=1.2, to=2, by=0.1)
        n_samp <- 1+rbinom(length(samp_times), 20, 0.2)
        mrt <- min(samp_times)
        Ne_times <- seq(from=1, to=20, length.out=1000)
        Ne_vals <- sapply(Ne_times, function(t) exp(0.3*t) * 0.1)
        phy_sim <- sim_coal(samp_times, n_samp, Ne_vals, Ne_times)
        tre <- build_coal_tree(phy_sim$samp_times, phy_sim$n_samp, phy_sim$coal_times)
        phy_ext <- extract_phy(tre)

        expect_equal(phy_sim$samp_times, phy_ext$samp_times+mrt)
        expect_equal(phy_sim$coal_times, phy_ext$coal_times+mrt)
        expect_equal(phy_sim$n_samp, phy_ext$n_samp)

        phy_offs <- offset_extracted(phy_ext, mrt)

        expect_equal(phy_sim$samp_times, phy_offs$samp_times)
        expect_equal(phy_sim$coal_times, phy_offs$coal_times)
        expect_equal(phy_sim$n_samp, phy_offs$n_samp)
    }
)

test_that("phy2poi correctly extracts coalescent times.",
{
    set.seed(34567)
    t_begin <- 0
    t_end <- 10

    samp_times <- seq(from=0, to=4, length.out=5)
    mrst <- min(samp_times)
    n_samp <- rnbinom(5, 5, 0.4)
    offset <- .1
    phy <- sim_coal(samp_times, n_samp, c(10,10), c(t_begin, t_end))

    tr <- build_coal_tree(phy$samp_times, phy$n_samp, phy$coal_times)

    ppp <- phy2poi(tr, mrst+offset, t_begin, t_end, max_dt=1e-2)

    At <-  compute_At(phy$samp_times,
                               phy$n_samp, phy$coal_times, t_end)

    coal_idx <- which(ppp$poi_counts > 0)
    process_times <- ppp$process_times

    process_times_scaled <- t_end-rescale_times(process_times, c(-1,1), c(t_begin, t_end))
    coal_times_scaled <- process_times_scaled[coal_idx]
    coal_times_scaled <- coal_times_scaled[order(coal_times_scaled)]
    phy_coal_truncated <- phy$coal_times[which(phy$coal_times <= t_end-t_begin)]+t_begin+offset
    
    expect_equal(coal_times_scaled, phy_coal_truncated)
})

test_that("phy2poi correctly agrees on precomputed phylogeny",
{
    t_begin <- 0
    t_end <- 4
    samp_times <- c(1,2,3.5)
    n_samp <- c(4,6,2)
    coal_times <- c(1.5, 2.5, 2.8, 3)

    At_vals <-   c(4,    4,   3,    3,   9,    9,   8,    8,    7,     6,    6,     8,    8)
    At_times <-  c(1.25, 1.5, 1.75, 2,   2.25, 2.5, 2.65, 2.8,  3,     3.25, 3.5,   3.75, 4)
    poi_emits <- c(0,    1,   0,    0,   0,    1,   0,    1,    1,     0,    0,     0,    0)
    wt        <- c(.25,  .25, .25,  .25, .25,  .25, .15,  .15,  .2,    .25,  .25,   .25,  .25)
    
    tr <- build_coal_tree(samp_times, n_samp, c(coal_times, seq(from=4.1, to=6, length.out=7)))
    ppp <- phy2poi(tr, 1.0, t_begin, t_end, max_dt=0.125)


    ppp$process_times <- t_end - rescale_times(ppp$process_times, c(-1,1), c(0,4))
    ppp$lineage_counts <- ppp$lineage_counts[order(ppp$process_times)]
    ppp$poi_counts <- ppp$poi_counts[order(ppp$process_times)]
    ppp$holding_times <- ppp$holding_times[order(ppp$process_times)]*2
    ppp$process_times  <- ppp$process_times [order(ppp$process_times)]

    expect_equal(ppp$process_times , At_times)
    expect_equal(ppp$lineage_counts, At_vals)
    expect_equal(ppp$poi_counts, poi_emits)
    expect_equal(ppp$holding_times, wt)
})

test_that("Approximate model likelihood equals exact likelihood for constant eff. pop. size",
{
    set.seed(34567)
    t_begin <- 0
    t_end <- 10

    span <- t_end-t_begin

    samp_times <- seq(from=0, to=4, length.out=5)
    mrst <- min(samp_times)
    n_samp <- rnbinom(5, 5, 0.4)
    phy <- sim_coal(samp_times, n_samp, c(10,10), c(t_begin, t_end))

    log_lh_coal <- 0
    inv_Ne <- 5
    s_idx <- 2
    c_idx <- 1
    n <- length(phy$samp_times)
    m <- length(phy$coal_times)
    t <- phy$samp_times[1]
    lins <- phy$n_samp[1]

    scale_factor <- (t_end-t_begin)/2

    while (c_idx <= m) {
        if (s_idx <= n && 
            phy$samp_times[s_idx] < phy$coal_times[c_idx]){

            t_next <- phy$samp_times[s_idx]
            lins_next <- lins + n_samp[s_idx]
            s_idx_next <- s_idx + 1
            c_idx_next <- c_idx
            sampling_evt <- TRUE
        } else {
            t_next <- phy$coal_times[c_idx]
            lins_next <- lins - 1
            s_idx_next <- s_idx
            c_idx_next <- c_idx + 1
            sampling_evt <- FALSE
        } 

        if (t_next < span) {
            log_lh_coal <- log_lh_coal - 0.5*lins*(lins-1)*(t_next-t)*inv_Ne
            if (!sampling_evt) {
                log_lh_coal <- log_lh_coal + log(0.5*lins*(lins-1)*inv_Ne)
            }
            s_idx <- s_idx_next
            c_idx <- c_idx_next
            lins <- lins_next
            t <- t_next
        } else {
            log_lh_coal <- log_lh_coal - 0.5*lins*(lins-1)*(span-t)*inv_Ne
            break
        }
    } 
    tr <- build_coal_tree(phy$samp_times, phy$n_samp, phy$coal_times)
    ppp <- phy2poi(tr, mrst, t_begin, t_end, max_dt=1e-3)

    n1 <- length(ppp$lineage_counts)

    rates_ppp <- sapply(c(1:n1), function(i) (log(0.5*(ppp$lineage_counts[i])*(ppp$lineage_counts[i]-1))+log(inv_Ne)))
    log_lh_ppp <- -sum(exp(rates_ppp+log(ppp$holding_times)+log(scale_factor))) + sum(rates_ppp[which(ppp$poi_counts > 0.5)])
    
    expect_equal(log_lh_ppp, log_lh_coal)
})

test_that("Approximate model likelihood is invariant to mesh refinement under constant rate",
{
    set.seed(34567)
    t_begin <- 0
    t_end <- 10

    samp_times <- seq(from=0, to=4, length.out=5)
    mrst <- min(samp_times)
    n_samp <- rnbinom(5, 5, 0.4)
    offset <- .1
    phy <- sim_coal(samp_times, n_samp, c(10,10), c(t_begin, t_end))

    tr <- build_coal_tree(phy$samp_times, phy$n_samp, phy$coal_times)

    ppp1 <- phy2poi(tr, mrst+offset, t_begin, t_end, max_dt=1e-3)
    ppp2 <- phy2poi(tr, mrst+offset, t_begin, t_end, max_dt=1e-4)

    n1 <- length(ppp1$lineage_counts)
    n2 <- length(ppp2$lineage_counts)
    C <- 5

    rates1 <- sapply(c(1:n1), function(i) log(0.5*(ppp1$lineage_counts[i])*(ppp1$lineage_counts[i]-1))+log(C))
    log_lh1 <- -sum(exp(rates1+log(ppp1$holding_times))) + sum(rates1[which(ppp1$poi_counts > 0.5)])

    rates2 <- sapply(c(1:n2), function(i) log(0.5*(ppp2$lineage_counts[i])*(ppp2$lineage_counts[i]-1))+log(C))
    log_lh2 <- -sum(exp(rates2+log(ppp2$holding_times))) + sum(rates2[which(ppp2$poi_counts > 0.5)])

    expect_equal(log_lh1, log_lh2)
})

test_that("Trajectory interpolation returns correct values", {
    usage_df <- simulated_usage()
    t0 <- min(usage_df$time)
    tmax <- max(usage_df$time)
    usage_f <- yearly_usg_stepfunc(usage_df$usage, usage_df$time)

    N <- 3e7 
    incidence <- 1.8e6 / N
    gamma_sus <- 1/12.0*365.0
    mu <- incidence+gamma_sus
    r0 <- mu/(gamma_sus)
    gamma_res_u <- 1/11.8*365.0
    gamma_res_t <- 1/12.3*365.0

    end_eq <- (1-1/r0)*N
    I_s0 <- 1*end_eq/11
    I_r0 <- 1*end_eq/11
    u0 <- list(S=N-I_s0-I_r0, I_s=I_s0, I_r=I_r0)

    beta_func <- function (t) if(t < (t0+(tmax-t0)/2)) mu*0.998 else mu * 1.002
    n_step <- 1e3
    df <- simulate_epi_nstrain(t0, tmax, n_step, n_strains=2, S_0 = unlist(u0[1]), I_0 = unlist(u0[2:3]), usage_f, beta_func, function(t) N, gamma_u=c(gamma_sus, gamma_res_u), gamma_t=c(gamma_sus, gamma_res_t))

    eval_times <- df$t 
    dt <- eval_times[2]-eval_times[1]

    span <- max(eval_times) - min(eval_times)

    k <- length(eval_times)
    Neg_s <- sapply((1:k), function(i) (df$I_s[i]^2)/(2*df$births_I_s[i]/dt))
    Neg_r <- sapply((1:k), function(i) (df$I_r[i]^2)/(2*df$births_I_r[i]/dt))

    test_times <- eval_times-min(eval_times)+0.5*dt ## times to test against.
                                                    ## Add 0.5dt so that these don't precisely correspond to knots 
                                                    ## Backwards time relative to present, i.e. 0 = present
                                                    ## Should have the f(test_times[i]) == Neg[k-i+1]

    f1 <- rev_time_step_func(Neg_s, t0, tmax, dt)
    f2 <- rev_time_step_func(Neg_r, t0, tmax, dt)

    a1 <- sapply((1:(k-1)), function(i) Neg_s[k-i+1])
    a2 <- sapply((1:(k-1)), function(i) Neg_r[k-i+1])

    s1 <- sapply(test_times[1:(k-1)], f1)
    s2 <- sapply(test_times[1:(k-1)], f2)

    expect_equal(a1, s1, tolerance=1e-10)
    expect_equal(a2, s2, tolerance=1e-10)
})

test_that("Simulation usage interpolation function agrees with model interpolation", {

    stan_mod <- cmdstan_model(system.file('stan',    
                                'test_interpolation.stan',
                                package='ResistPhy',
                                mustWork = T),
                            include_paths=system.file('stan',package='ResistPhy'),
                            compile=F) 

    stan_mod$compile()
    usage_df <- simulated_usage()

    N <- 1000
    t0 <- min(usage_df$time)
    tmax <- max(usage_df$time)

    t_eval <- seq(from=t0, to=tmax, length.out=N)

    usage_f <- yearly_usg_stepfunc(usage_df$usage, usage_df$time)

    dat <- list(M=length(usage_df$time),
                N=N,
                ABX_usage=usage_df$usage,
                usage_ts=usage_df$time,
                eval_ts=t_eval)

    stan_interp <- stan_mod$sample(data=dat,
            iter_warmup=0, iter_sampling=1, chains=1, seed=494838,
            fixed_param =T, refresh=1)

    ddf<-as_draws_df(stan_interp$draws())
    stan_res <- sapply(c(1:N), function(i) ddf[[paste0("usage_interp[",i,"]")]])
    r_res <- sapply(t_eval, usage_f)
    expect_equal(stan_res,r_res, tolerance=1e-6)
})