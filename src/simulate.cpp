#include <Rcpp.h>
#include <assert.h> 
using namespace Rcpp;

//' Simulate SIS epidemic process
// [[Rcpp::export]]
DataFrame simulate_fast(NumericVector u0,
            NumericVector usage_vals,
            NumericVector beta_vals,
            NumericVector times,
            double gamma_sus,
            double gamma_res_u,
            double gamma_res_t,
            double dt,
            double N) {

    int n = times.length();
    NumericVector S_v(n, 0.0);
    NumericVector I_s_v(n, 0.0);
    NumericVector I_r_v(n, 0.0);
    NumericVector births_I_s_v(n, 0.0);
    NumericVector births_I_r_v(n, 0.0);
    NumericVector inf_rate_v(n, 0.0);
    NumericVector rt_s_v(n,0.0);
    NumericVector rt_r_v(n,0.0);
    I_s_v[0] = u0[0];
    I_r_v[0] = u0[1];
    S_v[0] = u0[2];
    inf_rate_v[0] = beta_vals[0]*S_v[0]/N;

    rt_s_v[0] = inf_rate_v[0]/gamma_sus;
    rt_r_v[0] = inf_rate_v[0]/((1.0-usage_vals[0])*gamma_res_u + usage_vals[0]*gamma_res_t);
    
    for (int i = 1; i < n; i++) {
        double gamma_res = (1.0-usage_vals[i-1])*gamma_res_u + usage_vals[i-1]*gamma_res_t;
        double deaths_I_s = R::rpois(I_s_v[i-1]*gamma_sus*dt);
        double deaths_I_r = R::rpois(I_r_v[i-1]*gamma_res*dt);
        double births_I_s = R::rpois(I_s_v[i-1]*inf_rate_v[i-1]*dt);
        double births_I_r = R::rpois(I_r_v[i-1]*inf_rate_v[i-1]*dt);
        I_s_v[i] = I_s_v[i-1] + births_I_s-deaths_I_s;
        I_r_v[i] = I_r_v[i-1] + births_I_r-deaths_I_r;
        S_v[i] = S_v[i-1] + deaths_I_r + deaths_I_s - births_I_s - births_I_r;
        births_I_s_v[i] = births_I_s;
        births_I_r_v[i] = births_I_r;

        inf_rate_v[i] = beta_vals[i]*S_v[i]/N;
        rt_s_v[i] = inf_rate_v[i]/gamma_sus;
        rt_r_v[i] = inf_rate_v[i]/((1.0-usage_vals[i])*gamma_res_u + usage_vals[i]*gamma_res_t);
    }
    DataFrame out = DataFrame::create(Named("t") = times,         
                                      Named("I_s") = I_s_v,
                                      Named("I_r") = I_r_v,
                                      Named("S") = S_v,
                                      Named("births_I_s") = births_I_s_v,
                                      Named("births_I_r") = births_I_r_v,
                                      Named("inf_rate") = inf_rate_v,
                                      Named("Rt_s") = rt_s_v,
                                      Named("Rt_r") = rt_r_v);
    return  out;
}

//' Simulate SIS epidemic process
// [[Rcpp::export]]
DataFrame simulate_fast_nstrain(
            int n_strains,
            NumericVector I0,
            double S0,
            NumericVector usage_vals,
            NumericVector beta_vals,
            NumericVector N_vals,
            NumericVector gamma_u,
            NumericVector gamma_t,
            NumericVector times,
            double dt) {

    int n = times.length();

    NumericVector S_v(n, 0.0);
    NumericVector inf_rate_v(n, 0.0);
    std::vector<NumericVector> I_v;
    std::vector<NumericVector> births_v;
    std::vector<NumericVector> rt_v;

    S_v[0] = S0;
    inf_rate_v[0] = beta_vals[0]*S_v[0]/N_vals[0];

    for (int i = 0; i < n_strains; i++) {
        I_v.push_back(NumericVector(n, 0.0));
        I_v[i][0] = I0[i];
        births_v.push_back(NumericVector(n, 0.0));
        rt_v.push_back(NumericVector(n, 0.0));
        rt_v[i][0] = inf_rate_v[0]/((1.0-usage_vals[0])*gamma_u[i] + usage_vals[0]*gamma_t[i]);
    }
    
    for (int i = 1; i < n; i++) {
        std::vector<double> b_I(n_strains);
        std::vector<double> d_I(n_strains);
        std::vector<double> gamma(n_strains);

        for (int j = 0; j < n_strains; j++) {
            gamma[j] = ((1.0-usage_vals[i])*gamma_u[j] + usage_vals[i]*gamma_t[j]);
            b_I[j] = R::rpois(I_v[j][i-1]*inf_rate_v[i-1]*dt);
            d_I[j] = R::rpois(I_v[j][i-1]*gamma[j]*dt);
        }  

        S_v[i] = S_v[i-1] + std::accumulate(d_I.begin(), d_I.end(), 0.0) - std::accumulate(b_I.begin(), b_I.end(), 0.0); 
        S_v[i] += N_vals[i]-N_vals[i-1];
        inf_rate_v[i] = beta_vals[i]*S_v[i]/N_vals[i];

        for (int j = 0; j < n_strains; j++) {
            births_v[j][i] = b_I[j];
            I_v[j][i] = I_v[j][i-1] + b_I[j] - d_I[j];
            rt_v[j][i] = inf_rate_v[i]/gamma[j];
        }
    }

    CharacterVector namevec;
    DataFrame out;
    namevec.push_back("t");
    namevec.push_back("S");
    namevec.push_back("inf_rate");

    out.push_back(times);
    out.push_back(S_v);
    out.push_back(inf_rate_v);

    for (int i=0; i<n_strains; i++){
        int sidx = i+1;
        String name = String("I_");
        name+=String(sidx);
        namevec.push_back(name);
        out.push_back(I_v[i]);
    }     
    for (int i=0; i<n_strains; i++){
        int sidx = i+1;
        String name = String("births_");
        name+=String(sidx);
        namevec.push_back(name);
        out.push_back(births_v[i]);
    }    
    for (int i=0; i<n_strains; i++){
        int sidx = i+1;
        String name = String("Rt_");
        name+=String(sidx);
        namevec.push_back(name);
        out.push_back(rt_v[i]);
    }    

    out.attr("names") = namevec;
    Rcpp::DataFrame x;
    Rcpp::Language call("as.data.frame",out);
    x = call.eval();
    return x;
}


//' Simulate coalescent tree using thinning
// [[Rcpp::export]]
NumericVector sim_coal_fast(NumericVector samp_times,
            IntegerVector n_samp,
            NumericVector Ne_vals,
            NumericVector Ne_times,
            int samp_total,
            double dt,
            int n) {
    
    const double t_min = min(Ne_times);
    const double t_max = max(Ne_times);
    const double f_min = min(Ne_vals);
    const int m = samp_times.length();
    
    assert(Ne_times.length() == n);
    assert(Ne_vals.length() == n);
    assert(n_samp.length() == m);
    assert(dt > 0);

    NumericVector coal_times(samp_total-1, 0.0);

    std::function<double(double)> Ne_func = [Ne_vals, t_min, t_max, f_min, dt, n](double t) {
        assert(t >= 0.0);
        double out;
        double upper;
        upper = t_max-t_min;
        if (t >= upper) {
            out = f_min;
        } else {
            int index; 
            index = n-1-int(t/dt);
            assert(index < n);
            assert(index >= 0);
            out = Ne_vals[index];
        }
        return out;
    };

    int s_idx = 0;
    int c_idx = 0;
    int lin_count = n_samp[s_idx];
    double t_curr = samp_times[s_idx];

    s_idx++;

    while (lin_count > 1 || s_idx < m) {
        if (lin_count == 1) {
            t_curr = samp_times[s_idx];
            lin_count += n_samp[s_idx];
            s_idx++;
        } 
        assert(lin_count > 1);
        double s = R::rexp(f_min/(0.5*(lin_count)*(lin_count-1)));
        if (s_idx < m && t_curr+s >= samp_times[s_idx]) {
            t_curr = samp_times[s_idx];
            lin_count += n_samp[s_idx];
            s_idx++;
        } else {
            t_curr += s;
            double r = R::runif(0,1);
            if (r <= f_min/Ne_func(t_curr)) {
                coal_times[c_idx] = t_curr;
                lin_count--;
                c_idx++;
            }
        }
    }
    return coal_times;
}





