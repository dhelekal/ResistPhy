functions {
    #include basis_gp.stan
    #include utils.stan
}


data {
    int M; // Antibiotic usage array size
    int N_lineages; // Number of lineages. First Strain is always susceptible
    int obs_total; // Total number of observations across
    int coal_total;// total number of coalescent events 
    
    array[N_lineages] int<lower=0> obs_counts; // Observation counts for each lineage
    array[N_lineages] int<lower=0> coal_counts; // Coal event counts for each lineage
    
    vector[obs_total] times; // Times for all observations
    vector[obs_total] w_t; // Process holding times   

    array[coal_total] int<lower=0> coal_index; // Coalescent event index 
    array[obs_total] int<lower=0> lineage_counts; // A(t) values

    array[M] real usage_ts; // ABX usage time points
    array[M] real ABX_usage; // ABX usage values

    real<lower=0> time_scale; // GP Time Scale in days
    real<lower=0> unit_scale; // Time scale of the input 
    real<lower=0> gamma_guess; // Susceptible recovery rate
    real<lower=0> gamma_log_sd; // log sd for informative prior on gamma
    
    int K; // HSGP Approximation order
    real<lower=0> L; // HSGP Approximation boundary factor
}

transformed data {
    vector[obs_total] combNs = to_vector(choose(lineage_counts, 2));

    // Initialise  ragged basis matrix holding susceptible trajectory log-mean 
    matrix[obs_total,K] rate_accum_basis;
    for (k in 1:K) {
        for (j in 1:obs_total) {
            rate_accum_basis[j, k] = gp_basis_fun(times[j], L, k);
        }
    }

    // Initialise ragged basis matrix holding instantaneous rates for susceptible trajectory 
    matrix[obs_total,K] rate_inst_basis;
    for (k in 1:K) {
        for (j in 1:obs_total) {
            rate_inst_basis[j, k] = gp_basis_fun_deriv(times[j], L, k);
        }
    }

    /* Initialise ragged vectors holding constant function time integral,
     * usage function time integral, and usage function values. 
     * Initialise vectors of one timestep differences for these integrals.
     */
    vector[obs_total] const_int;
    vector[obs_total] usage_int;
    vector[obs_total] usage_vals;

    {
        int pos = 1;
        for (i in 1:N_lineages) {
            int N_obs = obs_counts[i];
            vector[N_obs] lineage_times = segment(times, pos, N_obs);

            vector[N_obs] temp_const_int = integrate_const(lineage_times);
            vector[N_obs] temp_usage_int = integrate_usage_step(lineage_times, ABX_usage, usage_ts, M);

            const_int[pos:(pos+N_obs-1)] = temp_const_int;
            usage_int[pos:(pos+N_obs-1)] = temp_usage_int;
            usage_vals[pos:(pos+N_obs-1)] = interp_usage_step(lineage_times, ABX_usage, usage_ts, N_obs, M);

            pos += N_obs;
        }
    }

    // Initialise scaling matrix used to decorrelate q-variables
    real q_scaling = 0.1*(1.0/4.0); //q parameter scaling to stabilise initialisation
    real usage_mean = mean(ABX_usage);
    matrix[2,2] q_rescale;
    q_rescale[1,:] = (1/(1-usage_mean))*[1.0, -1.0*usage_mean];
    q_rescale[2,:] = [0.0,1.0];
}

parameters {
    vector[K] f_tilde;
    vector[N_lineages] I_0_hat;
    matrix[N_lineages-1, 2] q_tilde;

    real gamma_sus_tilde; 
    real<lower=0> alpha;
    real<lower=0> rho; 
}

transformed parameters {
    real gamma_sus_sc = exp(gamma_sus_tilde*gamma_log_sd) * gamma_guess * time_scale;
    matrix[N_lineages-1, 2] q_log = q_scaling * q_tilde;
    vector[N_lineages-1] q_u;
    vector[N_lineages-1] q_t;

    vector[N_lineages] I_0_tilde = I_0_hat + 6;

    vector[K] sqrt_spd = diagSPD_EQ(alpha, rho, L, K);
    vector[K] coeffs = sqrt_spd .* f_tilde;

    vector[obs_total] traj_vec;
    vector[obs_total] log_mean_vec;
    vector[obs_total] log_beta_inst_vec;

    {
        array[N_lineages-1, 2] real q_arr = to_array_2d(exp(q_log) * q_rescale');

        vector[obs_total] beta_accum_vec = rate_accum_basis * coeffs;
        vector[obs_total] rate_inst_vec = rate_inst_basis * coeffs;

        log_beta_inst_vec = log(gamma_sus_sc) + log1p(rate_inst_vec/gamma_sus_sc);
        
        int pos_data = 1;
        for (i in 1:N_lineages) {
            int N_obs = obs_counts[i];
            vector[N_obs] lineage_times = segment(times, pos_data, N_obs);
            vector[N_obs] lineage_beta_accum = segment(beta_accum_vec, pos_data, N_obs);
            real c0 = I_0_tilde[i] - lineage_beta_accum[1];
            
            vector[N_obs] lineage_log_mean;

            if (i == 1) {
                lineage_log_mean = lineage_beta_accum + c0;
            } else {
                vector[N_obs] lineage_const_int = segment(const_int, pos_data, N_obs);                
                vector[N_obs] lineage_usage_int = segment(usage_int, pos_data, N_obs);
                q_u[i-1] = q_arr[i-1, 1];
                q_t[i-1] = q_arr[i-1, 2];

                lineage_log_mean = lineage_beta_accum + 
                                gamma_sus_sc * 
                                ((1-q_u[i-1]) * lineage_const_int - (q_t[i-1] - q_u[i-1]) * lineage_usage_int) + c0;
            }

            traj_vec[pos_data:(pos_data+N_obs-1)] = exp(lineage_log_mean);
            log_mean_vec[pos_data:(pos_data+N_obs-1)] = lineage_log_mean;
 
            pos_data += N_obs;
        }
    }
}
  
model {    
    f_tilde ~ normal(0, 1);
    gamma_sus_tilde ~ normal(0, 1); 
    q_u ~ lognormal(0,0.5);
    q_t ~ lognormal(0,0.5);
    alpha ~ gamma(4, 4);
    rho ~ inv_gamma(4.63,2.21);
    I_0_hat ~ normal(0,2);

    //Jacobian for q_tilde -> q_u, q_t transform

    target += sum(q_log);

    int pos_data = 1;
    int pos_coal = 1;
    for (i in 1:N_lineages) {

        int N_obs = obs_counts[i];

        vector[N_obs] lineage_log_traj = segment(log_mean_vec, pos_data, N_obs);
        vector[N_obs] lineage_log_beta_inst = segment(log_beta_inst_vec, pos_data, N_obs);
        vector[N_obs] lineage_combNs = segment(combNs, pos_data, N_obs);
        array[coal_counts[i]] int lineage_coal_index = segment(coal_index, pos_coal, coal_counts[i]);
        vector[N_obs] lineage_wt = segment(w_t, pos_data, N_obs);

        vector[N_obs] log_rates =  log(2.0) + lineage_log_beta_inst - lineage_log_traj;
        target += -sum(lineage_combNs .* lineage_wt .* exp(log_rates));
        target += sum(log_rates[lineage_coal_index]);

        pos_data += N_obs;
        pos_coal += coal_counts[i];
    }
}

generated quantities {
    real gamma_sus = gamma_sus_sc / time_scale;

    vector[N_lineages] I_0 = exp(I_0_tilde);
    vector[N_lineages-1] gamma_u = gamma_sus * q_u;
    vector[N_lineages-1] gamma_t = gamma_sus * q_t;
    
    vector[obs_total] birthrate = exp(log_beta_inst_vec)/time_scale * unit_scale;
    vector[obs_total] Ne = 0.5 * traj_vec ./ birthrate;
    vector[obs_total] I_mean = exp(log_mean_vec);
    vector[obs_total] R_t;

    {
        int pos = 1;
        for (i in 1:N_lineages) {
            int N_obs = obs_counts[i];
            vector[N_obs] lineage_birthrate = segment(birthrate, pos, N_obs);
            if (i==1) {
                R_t[pos:(pos+N_obs-1)] = lineage_birthrate ./ unit_scale ./ gamma_sus; 
            } else {
                vector[N_obs] lineage_usage_vals = segment(usage_vals, pos, N_obs);
                R_t[pos:(pos+N_obs-1)] = lineage_birthrate ./ 
                                        unit_scale ./ 
                                        (lineage_usage_vals*(gamma_t[i-1]-gamma_u[i-1]) + gamma_u[i-1]); 
            }
            pos += N_obs;
        }
    }
}

