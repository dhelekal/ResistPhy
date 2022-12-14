// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// simulate_fast
DataFrame simulate_fast(NumericVector u0, NumericVector usage_vals, NumericVector beta_vals, NumericVector times, double gamma_sus, double gamma_res_u, double gamma_res_t, double dt, double N);
RcppExport SEXP _ResistPhy_simulate_fast(SEXP u0SEXP, SEXP usage_valsSEXP, SEXP beta_valsSEXP, SEXP timesSEXP, SEXP gamma_susSEXP, SEXP gamma_res_uSEXP, SEXP gamma_res_tSEXP, SEXP dtSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type u0(u0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type usage_vals(usage_valsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_vals(beta_valsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_sus(gamma_susSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_res_u(gamma_res_uSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_res_t(gamma_res_tSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_fast(u0, usage_vals, beta_vals, times, gamma_sus, gamma_res_u, gamma_res_t, dt, N));
    return rcpp_result_gen;
END_RCPP
}
// simulate_fast_nstrain
DataFrame simulate_fast_nstrain(int n_strains, NumericVector I0, double S0, NumericVector usage_vals, NumericVector beta_vals, NumericVector N_vals, NumericVector gamma_u, NumericVector gamma_t, NumericVector times, double dt);
RcppExport SEXP _ResistPhy_simulate_fast_nstrain(SEXP n_strainsSEXP, SEXP I0SEXP, SEXP S0SEXP, SEXP usage_valsSEXP, SEXP beta_valsSEXP, SEXP N_valsSEXP, SEXP gamma_uSEXP, SEXP gamma_tSEXP, SEXP timesSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_strains(n_strainsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type I0(I0SEXP);
    Rcpp::traits::input_parameter< double >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type usage_vals(usage_valsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_vals(beta_valsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type N_vals(N_valsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma_u(gamma_uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma_t(gamma_tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_fast_nstrain(n_strains, I0, S0, usage_vals, beta_vals, N_vals, gamma_u, gamma_t, times, dt));
    return rcpp_result_gen;
END_RCPP
}
// sim_coal_fast
NumericVector sim_coal_fast(NumericVector samp_times, IntegerVector n_samp, NumericVector Ne_vals, NumericVector Ne_times, int samp_total, double dt, int n);
RcppExport SEXP _ResistPhy_sim_coal_fast(SEXP samp_timesSEXP, SEXP n_sampSEXP, SEXP Ne_valsSEXP, SEXP Ne_timesSEXP, SEXP samp_totalSEXP, SEXP dtSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type samp_times(samp_timesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_samp(n_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ne_vals(Ne_valsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ne_times(Ne_timesSEXP);
    Rcpp::traits::input_parameter< int >::type samp_total(samp_totalSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_coal_fast(samp_times, n_samp, Ne_vals, Ne_times, samp_total, dt, n));
    return rcpp_result_gen;
END_RCPP
}
// compute_At
DataFrame compute_At(NumericVector samp_times, IntegerVector n_samp, NumericVector coal_times, double t_max);
RcppExport SEXP _ResistPhy_compute_At(SEXP samp_timesSEXP, SEXP n_sampSEXP, SEXP coal_timesSEXP, SEXP t_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type samp_times(samp_timesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_samp(n_sampSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type coal_times(coal_timesSEXP);
    Rcpp::traits::input_parameter< double >::type t_max(t_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_At(samp_times, n_samp, coal_times, t_max));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ResistPhy_simulate_fast", (DL_FUNC) &_ResistPhy_simulate_fast, 9},
    {"_ResistPhy_simulate_fast_nstrain", (DL_FUNC) &_ResistPhy_simulate_fast_nstrain, 10},
    {"_ResistPhy_sim_coal_fast", (DL_FUNC) &_ResistPhy_sim_coal_fast, 7},
    {"_ResistPhy_compute_At", (DL_FUNC) &_ResistPhy_compute_At, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_ResistPhy(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
