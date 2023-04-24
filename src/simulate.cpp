#include <Rcpp.h>
#include <assert.h> 
#include <stdexcept>
#include <cstring>
using namespace Rcpp;

//' Simulate SIS epidemic process
//' @param u0 Vector u0
//' @param usage_vals Vector of usage values
//' @param beta_vals Vector of beta values
//' @param times Vector of times
//' @param gamma_sus Recovery rate for sus
//' @param gamma_res_u Recovery rate for res untreated
//' @param gamma_res_t Recovery rate for res treated
//' @param dt dt
//' @param N N
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
    I_s_v.at(0) = u0.at(0);
    I_r_v.at(0) = u0.at(1);
    S_v.at(0) = u0.at(2);
    inf_rate_v.at(0) = beta_vals.at(0)*S_v.at(0)/N;

    rt_s_v.at(0) = inf_rate_v.at(0)-gamma_sus;
    rt_r_v.at(0) = inf_rate_v.at(0)-((1.0-usage_vals.at(0))*gamma_res_u + usage_vals.at(0)*gamma_res_t);
    
    for (int i = 1; i < n; i++) {
        double gamma_res = (1.0-usage_vals.at(i-1))*gamma_res_u + usage_vals.at(i-1)*gamma_res_t;
        double deaths_I_s = R::rpois(I_s_v.at(i-1)*gamma_sus*dt);
        double deaths_I_r = R::rpois(I_r_v.at(i-1)*gamma_res*dt);
        double births_I_s = R::rpois(I_s_v.at(i-1)*inf_rate_v.at(i-1)*dt);
        double births_I_r = R::rpois(I_r_v.at(i-1)*inf_rate_v.at(i-1)*dt);
        I_s_v.at(i) = I_s_v.at(i-1) + births_I_s-deaths_I_s;
        I_r_v.at(i) = I_r_v.at(i-1) + births_I_r-deaths_I_r;
        S_v.at(i) = S_v.at(i-1) + deaths_I_r + deaths_I_s - births_I_s - births_I_r;
        births_I_s_v.at(i) = births_I_s;
        births_I_r_v.at(i) = births_I_r;

        inf_rate_v.at(i) = beta_vals.at(i)*S_v.at(i)/N;
        rt_s_v.at(i) = inf_rate_v.at(i)-gamma_sus;
        rt_r_v.at(i) = inf_rate_v.at(i)-((1.0-usage_vals.at(i))*gamma_res_u + usage_vals.at(i)*gamma_res_t);
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
//' @param n_strains Number of strains
//' @param I0 Vector of initial infecteds
//' @param S0 Initial susceptibles
//' @param usage_vals Vector of usage values
//' @param beta_vals Vector of beta values
//' @param N_vals Vector N_vals
//' @param gamma_u recovery rates when not treated with ABX of interest
//' @param gamma_t recovery rates when treated with ABX of interest
//' @param times Vector of times
//' @param dt dt
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

    S_v.at(0) = S0;
    inf_rate_v.at(0) = beta_vals.at(0)*S_v.at(0)/N_vals.at(0);

    for (int i = 0; i < n_strains; i++) {
        I_v.push_back(NumericVector(n, 0.0));
        I_v.at(i).at(0) = I0.at(i);
        births_v.push_back(NumericVector(n, 0.0));
        rt_v.push_back(NumericVector(n, 0.0));
        rt_v.at(i).at(0) = inf_rate_v.at(0)-((1.0-usage_vals.at(0))*gamma_u.at(i) + usage_vals.at(0)*gamma_t.at(i));
    }
    
    for (int i = 1; i < n; i++) {
        std::vector<double> b_I(n_strains);
        std::vector<double> d_I(n_strains);
        std::vector<double> gamma(n_strains);

        for (int j = 0; j < n_strains; j++) {
            gamma.at(j) = ((1.0-usage_vals.at(i))*gamma_u.at(j) + usage_vals.at(i)*gamma_t.at(j));
            b_I.at(j) = R::rpois(I_v.at(j).at(i-1)*inf_rate_v.at(i-1)*dt);
            d_I.at(j) = R::rpois(I_v.at(j).at(i-1)*gamma.at(j)*dt);
        }  

        S_v.at(i) = S_v.at(i-1) + std::accumulate(d_I.begin(), d_I.end(), 0.0) - std::accumulate(b_I.begin(), b_I.end(), 0.0); 
        S_v.at(i) += N_vals.at(i)-N_vals.at(i-1);
        inf_rate_v.at(i) = beta_vals.at(i)*S_v.at(i)/N_vals.at(i);

        for (int j = 0; j < n_strains; j++) {
            births_v.at(j).at(i) = b_I.at(j);
            I_v.at(j).at(i) = I_v.at(j).at(i-1) + b_I.at(j) - d_I.at(j);
            rt_v.at(j).at(i) = inf_rate_v.at(i)-gamma.at(j);
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
        out.push_back(I_v.at(i));
    }     
    for (int i=0; i<n_strains; i++){
        int sidx = i+1;
        String name = String("births_");
        name+=String(sidx);
        namevec.push_back(name);
        out.push_back(births_v.at(i));
    }    
    for (int i=0; i<n_strains; i++){
        int sidx = i+1;
        String name = String("Rt_");
        name+=String(sidx);
        namevec.push_back(name);
        out.push_back(rt_v.at(i));
    }    

    out.attr("names") = namevec;
    Rcpp::DataFrame x;
    Rcpp::Language call("as.data.frame",out);
    x = call.eval();
    return x;
}


//' Simulate coalescent tree using thinning
//' @param samp_times a vector of sampling times
//' @param n_samp a vector containing the numbers of samples taken at corresponding sampling time.
//' @param Ne_vals a vector of Ne(t) values.
//' @param Ne_times a vector of time points corresponding to Ne(t) values.
//' @param samp_total samp_total
//' @param dt dt
//' @param n n
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

/*    std::function<double(double)> Ne_func = [Ne_vals, t_min, t_max, f_min, dt, n](double t) {
        assert(t >= 0.0);
        double out;
        const double upper = t_max-t_min;
        if (t >= upper) {
            out = f_min;
        } else {
            const int index = n-1-int(t/dt);
            if (index >= n || index < 0) throw 
                std::logic_error("Ne index out of bounds! Index: " + std::to_string(index) + " Extent: " + std::to_string(n) + " t: " +std::to_string(t) + " dt: " + std::to_string(dt)); 
            assert(index < n);
            assert(index >= 0);
            out = Ne_vals.at(index);
        }
        return out;
    };*/

    int ne_idx = n-2;

    int s_idx = 0;
    int c_idx = 0;
    int lin_count = n_samp.at(s_idx);
    double t_curr = samp_times.at(s_idx);
    
    const double upper = t_max-t_min;

    s_idx++;

    while (lin_count > 1 || s_idx < m) {
        if (lin_count == 1) {
            t_curr = samp_times.at(s_idx);
            lin_count += n_samp.at(s_idx);
            s_idx++;
        } 
        assert(lin_count > 1);
        double s = R::rexp(f_min/(0.5*(lin_count)*(lin_count-1)));
        if (s_idx < m && t_curr+s >= samp_times.at(s_idx)) {
            t_curr = samp_times.at(s_idx);
            lin_count += n_samp.at(s_idx);
            s_idx++;
        } else {
            t_curr += s;
            
            double ne_val;

            if (t_curr >= upper)
            {
                ne_val = f_min;
            } else 
            {
                const double tr = t_max-t_curr;
                while (ne_idx > 1 && tr < Ne_times.at(ne_idx))
                {
                    ne_idx--;
                }
                
                if (!(tr > Ne_times.at(ne_idx) && tr <= Ne_times.at(ne_idx+1)))
                    std::logic_error("Ne index localisation failed: " + std::to_string(ne_idx)); 
                
                ne_val = Ne_vals.at(ne_idx+1);
            }

            double r = R::runif(0.0,1.0);
            if (r <= f_min/ne_val) {
                coal_times.at(c_idx) = t_curr;
                lin_count--;
                c_idx++;
            }
        }
    }
    return coal_times;
}





