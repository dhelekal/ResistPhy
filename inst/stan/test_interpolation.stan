functions {
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
    int M;
    int N;
    array[M] real ABX_usage; // Antibiotic Usage values
    array[M] real usage_ts; // Antibiotic Usage time points
    vector[N] eval_ts;
}

generated quantities {
    vector[N] usage_interp = interp_usage_step(eval_ts, ABX_usage, usage_ts, N, M);
}