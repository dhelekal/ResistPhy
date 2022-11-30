
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

    real gp_basis_fun_int(real t, real s, real L, int k) {
        real a = -2*sqrt(L)/(k*pi());
        real hi = cos(sqrt_gp_lambda(k,L)*(L+s));
        real lo = cos(sqrt_gp_lambda(k,L)*(L+t));
        real res = a*(hi-lo);
        return res;
    }

    real spd_matern(real x, real alpha, real rho) {
       real dens = 4 * (alpha^2) * (sqrt(3)/rho)^3 * 1/((sqrt(3)/rho)^2 + x^2)^2;
       return dens;
    }

    real sqrt_spd_matern(real x, real alpha, real rho) {
       real dens = 2 * alpha * ((sqrt(3)/rho)^1.5) * 1/((sqrt(3)/rho)^2 + x^2);
       return dens;
    }

    vector vec_sqrt_spd_matern32(real alpha, real rho, real L, int K) {
         return 2 * alpha * ((sqrt(3)/rho)^1.5) * inv((sqrt(3)/rho)^2 + ((pi()/(2*L)) * linspaced_vector(K, 1, K))^2); //Adapted from aki vektharis case study
    }

    //vector vec_sqrt_spd_matern52(real rho, real L, int K) {
    //     return sqrt(16.0/3.0) * ((sqrt(5.0)/rho)^2.5) * inv(((sqrt(5.0)/rho)^2.0 + ((pi()/(2.0*L)) * linspaced_vector(K, 1, K))^2.0)^(1.5)); //Adapted from aki vektharis case study
    //}

    vector diagSPD_EQ(real alpha, real rho, real L, int M) {
        return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
    }

    vector vec_sqrt_spd_matern52(real alpha, real rho, real L, int K) {
        return 4 * sqrt(rho*((5.0^2.5)/3.0)) * alpha * inv((5.0 + ((rho*pi()/(2*L)) * linspaced_vector(K, 1, K))^2.0)^1.5); //Adapted from aki vektharis case study
    }

    //real vec_sqrt_spd_RBF(real rho,  real L, int K) {
    //		return sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(K, 1, K)^2);
    //}