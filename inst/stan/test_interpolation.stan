functions {
    #include utils.stan
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