functions {
    vector qdiff(vector y, vector theta, real[] x_r, int[] x_i) {
        vector[2] diff;
        diff[1] = inv_gamma_cdf(theta[1] | exp(y[1]), exp(y[2])) - 0.01;
        diff[2] = 1 - inv_gamma_cdf(theta[2] | exp(y[1]), exp(y[2])) - 0.01;
        return diff;
    }
}

transformed data {
    real lo = 0.2; 
    real hi = 2;   

    array[0] real x_r;
    array[0] int x_i;
    vector[2] y = algebra_solver(qdiff, [1,1]', [lo,hi]', x_r, x_i);
}

generated quantities {
    real alpha = exp(y[1]);
    real beta = exp(y[2]);
}