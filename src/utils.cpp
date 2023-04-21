#include <Rcpp.h>
#include <assert.h> 
using namespace Rcpp;

//' Compute A(t), lineages through time
//' @param samp_times a vector of sampling times
//' @param n_samp a vector containing the numbers of samples taken at corresponding sampling time
//' @param coal_times vector of coalescent times
//' @param t_max maximum time
// [[Rcpp::export]]

DataFrame compute_At(NumericVector samp_times,
                    IntegerVector n_samp,
                    NumericVector coal_times,
                    double t_max){
    DataFrame out;

    NumericVector vals = NumericVector::create();
    NumericVector times = NumericVector::create();

    bool in_bounds = true;
    long int s_idx = 0;
    long int c_idx = 0;
    double lin_count = 0; 

    while (in_bounds && (s_idx < samp_times.length() || c_idx < coal_times.length())) {
        double next_t;
        double prev_val = lin_count;
        if ((s_idx < samp_times.length()) && 
            (c_idx < coal_times.length()) && 
            (samp_times.at(s_idx) <= coal_times.at(c_idx))) {
            
            next_t = samp_times.at(s_idx);
            lin_count += n_samp.at(s_idx);
            s_idx++;
        } else if ((s_idx < samp_times.length()) && 
                (c_idx < coal_times.length()) && 
                (samp_times.at(s_idx) > coal_times.at(c_idx))) {
            
            next_t = coal_times.at(c_idx);
            lin_count -= 1;
            c_idx++;
        } else if (s_idx < samp_times.length()) {

            next_t = samp_times.at(s_idx);
            lin_count += n_samp.at(s_idx);
            s_idx++;
        } else {

            next_t = coal_times.at(c_idx);
            lin_count -= 1;
            c_idx++; 
        }
        
        if (next_t <= t_max) {
            times.push_back(next_t);
            vals.push_back(lin_count);
        } else {
            times.push_back(t_max);
            vals.push_back(prev_val);
            in_bounds=false;
        }
    }

    out = DataFrame::create(Named("t") = times,         
                        Named("At") = vals);
    return out;
}