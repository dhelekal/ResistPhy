    vector compute_combNs(array[] int lin_count, int N) {
        vector[N] out;  
        for (i in 1:N) {
            out[i] = 0.5 * (lin_count[i]) * (lin_count[i]-1);
        }
        return(out);
    }

    real tree_log_lh(int K, 
                    array[] int event_type,
                    vector accum_rate,
                    vector log_inv_Neg,
                    vector combNs,
                    real scale) {
        vector[K] summands;
        summands[1] = 0.0; //First event is always sampling with no lineages present
        for (i in 2:K) {
            summands[i] = -combNs[i] * (accum_rate[i] - accum_rate[i-1])*scale;
            if (event_type[i] == 1) {
                if (is_inf(log(combNs[i]))) reject("log combn is inf");
                summands[i] += log_inv_Neg[i] + log(scale) + log(combNs[i]);
            }
        }
        return(sum(summands));
    }

    real tree_log_lh2(int K, 
                    array[] int event_type,
                    vector interval_rates,
                    vector log_inv_Neg,
                    vector combNs,
                    real scale) {
        vector[K] summands;
        summands[1] = 0.0; //First event is always sampling with no lineages present
        for (i in 2:K) {
            if (interval_rates[i] < 0) reject("Invalid rate!"); 
            summands[i] = -combNs[i] * interval_rates[i] * scale;
            if (event_type[i] == 1) {
                if (is_inf(log(combNs[i]))) reject("Log combn is inf!");
                summands[i] += log_inv_Neg[i] + log(scale) + log(combNs[i]);
            }
        }
        return(sum(summands));
    }    