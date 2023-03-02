functions {

    real sqrt_gp_lambda(int k, real L) {
        return (k*pi())/(2*L);
    }

    real gp_basis_fun(real t, real L, int k) {
        real S = 1/sqrt(L) * sin(sqrt_gp_lambda(k,L)*(t+L));
        return S;
    }

    real gp_basis_fun_deriv(real t, real L, int k) {
        return 1/sqrt(L) * sqrt_gp_lambda(k,L) * cos(sqrt_gp_lambda(k,L)*(t+L));
    }
    
    vector diagSPD_EQ(real alpha, real rho, real L, int M) {
        return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
    }

    real const_primitive(real t, real t_lower, real val) {
        real s = t - t_lower;
        real out = s*val;
        return out;
    }
    
    real negexp_primitive(real t, real t_lower, real rate) {
        real t_min = -1.0;
        real out = 1/rate * (exp(-rate * (t_lower - t_min)) - exp(-rate * (t - t_min)));
        return out;
    }

    //integrate linear interpolation of usage function, i.e. \int_{-1}^{t_i}{ u(t) }\,dt 
    vector integrate_usage_step(vector times, array[] real ABX_usage, array[] real usage_ts, int M){
        
        //First compute the integrals of all components
        vector[M-1] component_ints = rep_vector(0.0, M-1); 
        for(j in 2:M) {
            real t_begin = usage_ts[j-1];
            real usage_curr = ABX_usage[j-1];
            real t_end = usage_ts[j];

            component_ints[j-1] = const_primitive(t_end, t_begin, usage_curr);
        }
        
        int N = size(times);
        vector[N] out = rep_vector(0.0, N);

        for (i in 2:N) {
            real t = times[i];
            for(j in 2:M) {
                if (t > usage_ts[j-1] && t <= usage_ts[j]){
                    out[i] += const_primitive(t, usage_ts[j-1], ABX_usage[j-1]);
                    break;
                } else if (t > usage_ts[j-1]){
                    out[i] += component_ints[j-1];                    
                } else {
                    print(" ", t, ";");
                    reject("Illegal state reached");   
                }
            }
        }
        return out;
    }

    vector integrate_usage_decay(vector times, array[] real ABX_usage, array[] real usage_ts, int M) {
                
        //First compute the integrals of all components
        vector[M-1] component_ints = rep_vector(0.0, M-1);
        for(j in 2:M) {
            real t_begin = usage_ts[j-1];
            real usage_curr = ABX_usage[j-1];
            real t_end = usage_ts[j];

            component_ints[j-1] = const_primitive(t_end, t_begin, usage_curr);
        }
        
        int N = size(times);
        vector[N] out = rep_vector(0.0, N);

        for (i in 2:N) {
            real t = times[i];
            for(j in 2:M) {
                if (t > usage_ts[j-1] && t <= usage_ts[j]){
                    out[i] += const_primitive(t, usage_ts[j-1], ABX_usage[j-1]);
                    break;
                } else if (t > usage_ts[j-1]){
                    out[i] += component_ints[j-1];                    
                } else {
                    print(" ", t, ";");
                    reject("Illegal state reached");   
                }
            }
        }
        return out;
    }

    vector integrate_components_decay_usage(array[] real ABX_usage, array[] real usage_ts, real rate){
        //First compute the integrals of all components
        int M = size(usage_ts);
        vector[M-1] component_ints = rep_vector(0.0, M-1);
        for(j in 2:M) {
            real t_begin = usage_ts[j-1];
            real usage_curr = ABX_usage[j-1];
            real t_end = usage_ts[j];

            component_ints[j-1] = usage_curr * negexp_primitive(t_end, t_begin, rate);
        }
        return component_ints;
    } 

    vector integrate_components_decay_const(array[] real usage_ts, real rate){
        //First compute the integrals of all components
        int M = size(usage_ts);
        vector[M-1] component_ints = rep_vector(0.0, M-1);
        for(j in 2:M) {
            real t_begin = usage_ts[j-1];
            real t_end = usage_ts[j];

            component_ints[j-1] = negexp_primitive(t_end, t_begin, rate);
        }
        return component_ints;
    }

    vector localise_times_usage(vector times, array[] real ABX_usage, array[] real usage_ts){
        int N = size(times);
        int M = size(usage_ts);
        vector[N] out = rep_vector(0.0, N);
        for (i in 2:N) {
            real t = times[i];
            for(j in 2:M) {
                if (t > usage_ts[j-1] && t <= usage_ts[j]){
                    out[i] = ABX_usage[j-1];
                    break;
                } 
            }
        }
        return out;
    } 

    vector localise_times_lower(vector times, array[] real usage_ts){
        int N = size(times);
        int M = size(usage_ts);
        vector[N] out = rep_vector(0.0, N);
        for (i in 2:N) {
            real t = times[i];
            for(j in 2:M) {
                if (t > usage_ts[j-1] && t <= usage_ts[j]){
                    out[i] = usage_ts[j-1];
                    break;
                } 
            }
        }
        return out;
    } 

    
    matrix build_component_sum_mat(vector times, array[] real usage_ts){
        int N = size(times);
        int M = size(usage_ts);
        matrix[N,M-1] out = rep_matrix(0.0,N,M-1);
        for (i in 2:N) {
            real t = times[i];
            row_vector[M-1] r = rep_row_vector(0.0,M-1);
            for(j in 2:M) {
                if (t > usage_ts[j]){
                    r[j-1] = 1.0;
                } 
            }
            out[i, :] = r;
        }
        return out;
    } 

    //integrate constant function
    vector integrate_const(vector times) {
        int N = size(times);
        vector[N] out;

        for(i in 1:N) {
            out[i] = times[i]+1.0;
        }
        return out;
    }

    vector interp_usage_step(vector times, array[] real ABX_usage, array[] real usage_ts, int N, int M) {
        vector[N] out;
        for (i in 1:N) {
            real t = times[i];
            for(j in 1:(M-1)) {
                if (t >= usage_ts[j] && t < usage_ts[j+1]){
                    out[i] = ABX_usage[j];
                } else if (t >= usage_ts[M]) {
                    out[i] = ABX_usage[M];
                }
            }
        }
        return out;
    }

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

    int sus_count = obs_counts[1];
    int res_count = obs_total - sus_count;
    vector[res_count] res_Ts = segment(times, 1+sus_count, sum(obs_counts[2:N_lineages]));
    //Initialise usage time localisation for cost decay for resistant strains
    vector[res_count] time_loc = localise_times_lower(res_Ts, usage_ts);
    //Initialise usage value localisation for cost decay for resistant strains
    vector[res_count] usage_loc = localise_times_usage(res_Ts, ABX_usage, usage_ts);
    //build decay rate integral sum matrix for resistant strains
    matrix[res_count, M-1] comp_sum_mat = build_component_sum_mat(res_Ts, usage_ts);
    
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
    real usage_mean = mean(ABX_usage);
    matrix[2,2] q_rescale;
    real A = (1.0/(1.0-usage_mean));
    q_rescale[1,:] = [1.0,0.0];//A*[1.0, -1.0*usage_mean];
    q_rescale[2,:] = [0.0,1.0];
}

parameters {
    vector[K] f_tilde;
    vector[N_lineages] I_0_hat;
    matrix[N_lineages-1, 2] q_tilde;

    real gamma_sus_tilde; 
    real<lower=0> alpha;
    real<lower=0> rho; 

    array[N_lineages-1] real<lower=0> phi;
    array[N_lineages-1] real<lower =0, upper= 1> p;
}

transformed parameters {
    real gamma_sus_sc = exp(gamma_sus_tilde*gamma_log_sd) * gamma_guess * time_scale;
    matrix[N_lineages-1, 2] q_hat  = q_tilde * q_rescale' + gamma_sus_sc; //Is this how the rescaling works? verify.
    array[N_lineages-1, 2] real gamma_res_q = to_array_2d(log1p_exp(q_hat));

    vector[N_lineages-1] gamma_u_sc;
    vector[N_lineages-1] gamma_t_sc;

    vector[N_lineages] I_0_tilde = I_0_hat + 6;

    vector[K] sqrt_spd = diagSPD_EQ(alpha, rho, L, K);
    vector[K] coeffs = sqrt_spd .* f_tilde;

    vector[obs_total] traj_vec;
    vector[obs_total] log_mean_vec;
    
    vector[obs_total] log_beta_inst_vec;

    vector[res_count] log_mean_vec_decay;
    vector[res_count] traj_vec_decay;

    {
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
                vector[N_obs] lineage_log_mean_decay;

                int pos_res = pos_data-sus_count;
                vector[M-1] usage_decay_steps = integrate_components_decay_usage(ABX_usage, usage_ts, phi[i-1]);
                vector[M-1] const_decay_steps = integrate_components_decay_const(usage_ts, phi[i-1]);

                vector[N_obs] usage_decay_ints;
                vector[N_obs] const_decay_ints; 

                vector[N_obs] lineage_usage_loc = segment(usage_loc, pos_res, N_obs);
                vector[N_obs] lineage_time_loc = segment(time_loc, pos_res, N_obs);
                matrix[N_obs, M-1] lineage_sum_mat = block(comp_sum_mat, pos_res, 1,  N_obs, M-1);
        
                for (j in 1:N_obs) {
                    const_decay_ints[j] = negexp_primitive(lineage_times[j], lineage_time_loc[j], phi[i-1]);
                }

                for (j in 1:N_obs) {
                    usage_decay_ints[j] = lineage_usage_loc[i]*negexp_primitive(lineage_times[j], lineage_time_loc[j], phi[i-1]);
                }

                const_decay_ints += lineage_sum_mat * const_decay_steps;
                usage_decay_ints += lineage_sum_mat * usage_decay_steps;

                vector[N_obs] lineage_const_int = segment(const_int, pos_data, N_obs);                
                vector[N_obs] lineage_usage_int = segment(usage_int, pos_data, N_obs);

                gamma_u_sc[i-1] = gamma_res_q[i-1, 1];
                gamma_t_sc[i-1] = gamma_res_q[i-1, 2];

                lineage_log_mean = lineage_beta_accum + 
                                (gamma_sus_sc - gamma_u_sc[i-1]) * lineage_const_int - (gamma_t_sc[i-1] - gamma_u_sc[i-1]) * lineage_usage_int + c0;
                
                lineage_log_mean_decay = lineage_beta_accum + gamma_sus_sc * lineage_const_int -
                    gamma_t_sc[i-1] * lineage_usage_int -
                    (gamma_u_sc[i-1]-gamma_sus_sc) * (const_decay_ints - usage_decay_ints) - 
                    (lineage_const_int - lineage_usage_int) * gamma_sus_sc +
                    c0;

                traj_vec_decay[pos_res:(pos_res+N_obs-1)] = exp(lineage_log_mean_decay);
                log_mean_vec_decay[pos_res:(pos_res+N_obs-1)] = lineage_log_mean_decay;
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
    to_array_1d(gamma_res_q) ~ normal(gamma_sus_sc, 0.3 * gamma_guess * time_scale);
    alpha ~ gamma(4, 4);
    rho ~ inv_gamma(4.63,2.21);
    I_0_hat ~ normal(0,2);
    phi ~ normal(0, 1);

    //Jacobian for q_tilde -> gamma_res_q transform
    target += sum(log_inv_logit(q_hat));

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

        if (i == 1) {
            target += -sum(lineage_combNs .* lineage_wt .* exp(log_rates));
            target += sum(log_rates[lineage_coal_index]);
        } else {
            int pos_res = pos_data-sus_count;
            vector[N_obs] lineage_log_traj_decay = segment(log_mean_vec_decay, pos_res, N_obs);
            vector[N_obs] log_rates_decay = log(2.0) + lineage_log_beta_inst - lineage_log_traj_decay;
            
            real nodecay_inc = sum(log_rates[lineage_coal_index]) - sum(lineage_combNs .* lineage_wt .* exp(log_rates));
            real decay_inc = sum(log_rates_decay[lineage_coal_index]) - sum(lineage_combNs .* lineage_wt .* exp(log_rates_decay));

            target += log_sum_exp(log(p[i-1]) +  nodecay_inc, log1m(p[i-1]) + decay_inc);
        }

        pos_data += N_obs;
        pos_coal += coal_counts[i];
    }
}

generated quantities {
    real gamma_sus = gamma_sus_sc / time_scale * unit_scale;

    vector[N_lineages] I_0 = exp(I_0_tilde);
    vector[N_lineages-1] gamma_u = gamma_u_sc / time_scale * unit_scale;
    vector[N_lineages-1] gamma_t = gamma_t_sc / time_scale * unit_scale;
    vector[N_lineages-1] q_u = gamma_u - gamma_sus;
    vector[N_lineages-1] q_t = gamma_t - gamma_sus;
    
    vector[obs_total] birthrate = exp(log_beta_inst_vec)/time_scale * unit_scale;
    vector[obs_total] Ne = 0.5 * traj_vec ./ birthrate;
    vector[obs_total] I_mean = exp(log_mean_vec);
    vector[obs_total] r_t;

    {
        int pos = 1;
        for (i in 1:N_lineages) {
            int N_obs = obs_counts[i];
            vector[N_obs] lineage_birthrate = segment(birthrate, pos, N_obs);
            if (i==1) {
                r_t[pos:(pos+N_obs-1)] = lineage_birthrate ./ unit_scale  -  gamma_sus; 
            } else {
                vector[N_obs] lineage_usage_vals = segment(usage_vals, pos, N_obs);
                r_t[pos:(pos+N_obs-1)] = lineage_birthrate ./ 
                                        unit_scale - (lineage_usage_vals*(gamma_t[i-1]-gamma_u[i-1]) + gamma_u[i-1]); 
            }
            pos += N_obs;
        }
    }
}

