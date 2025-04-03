// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Kmatrix_
List Kmatrix_(int N, double Pr1, double Pr2, double Pr3, double Pr4);
RcppExport SEXP _SMCABCnJRNMM_Kmatrix_(SEXP NSEXP, SEXP Pr1SEXP, SEXP Pr2SEXP, SEXP Pr3SEXP, SEXP Pr4SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type Pr1(Pr1SEXP);
    Rcpp::traits::input_parameter< double >::type Pr2(Pr2SEXP);
    Rcpp::traits::input_parameter< double >::type Pr3(Pr3SEXP);
    Rcpp::traits::input_parameter< double >::type Pr4(Pr4SEXP);
    rcpp_result_gen = Rcpp::wrap(Kmatrix_(N, Pr1, Pr2, Pr3, Pr4));
    return rcpp_result_gen;
END_RCPP
}
// KmatrixgivenLc_
NumericMatrix KmatrixgivenLc_(int N, double L, double c);
RcppExport SEXP _SMCABCnJRNMM_KmatrixgivenLc_(SEXP NSEXP, SEXP LSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(KmatrixgivenLc_(N, L, c));
    return rcpp_result_gen;
END_RCPP
}
// perturb_continuous_withinprior_
NumericVector perturb_continuous_withinprior_(NumericVector theta_c_sampled, NumericMatrix sigma_kernel, NumericMatrix Pr_cont);
RcppExport SEXP _SMCABCnJRNMM_perturb_continuous_withinprior_(SEXP theta_c_sampledSEXP, SEXP sigma_kernelSEXP, SEXP Pr_contSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta_c_sampled(theta_c_sampledSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma_kernel(sigma_kernelSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pr_cont(Pr_contSEXP);
    rcpp_result_gen = Rcpp::wrap(perturb_continuous_withinprior_(theta_c_sampled, sigma_kernel, Pr_cont));
    return rcpp_result_gen;
END_RCPP
}
// perturb_continuous_
NumericVector perturb_continuous_(NumericVector theta_c_sampled, NumericMatrix sigma_kernel);
RcppExport SEXP _SMCABCnJRNMM_perturb_continuous_(SEXP theta_c_sampledSEXP, SEXP sigma_kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta_c_sampled(theta_c_sampledSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma_kernel(sigma_kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(perturb_continuous_(theta_c_sampled, sigma_kernel));
    return rcpp_result_gen;
END_RCPP
}
// nJRNMM_prior_
NumericVector nJRNMM_prior_(NumericVector theta, int draw, NumericMatrix Pr_cont);
RcppExport SEXP _SMCABCnJRNMM_nJRNMM_prior_(SEXP thetaSEXP, SEXP drawSEXP, SEXP Pr_contSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type draw(drawSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pr_cont(Pr_contSEXP);
    rcpp_result_gen = Rcpp::wrap(nJRNMM_prior_(theta, draw, Pr_cont));
    return rcpp_result_gen;
END_RCPP
}
// model_
NumericMatrix model_(int N, NumericMatrix sol);
RcppExport SEXP _SMCABCnJRNMM_model_(SEXP NSEXP, SEXP solSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sol(solSEXP);
    rcpp_result_gen = Rcpp::wrap(model_(N, sol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SMCABCnJRNMM_Kmatrix_", (DL_FUNC) &_SMCABCnJRNMM_Kmatrix_, 5},
    {"_SMCABCnJRNMM_KmatrixgivenLc_", (DL_FUNC) &_SMCABCnJRNMM_KmatrixgivenLc_, 3},
    {"_SMCABCnJRNMM_perturb_continuous_withinprior_", (DL_FUNC) &_SMCABCnJRNMM_perturb_continuous_withinprior_, 3},
    {"_SMCABCnJRNMM_perturb_continuous_", (DL_FUNC) &_SMCABCnJRNMM_perturb_continuous_, 2},
    {"_SMCABCnJRNMM_nJRNMM_prior_", (DL_FUNC) &_SMCABCnJRNMM_nJRNMM_prior_, 3},
    {"_SMCABCnJRNMM_model_", (DL_FUNC) &_SMCABCnJRNMM_model_, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SMCABCnJRNMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
