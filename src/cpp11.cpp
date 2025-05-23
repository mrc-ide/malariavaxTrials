// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// acquire_immunity.cpp
doubles acquire_immunity_cpp(doubles exposure, double u, double d);
extern "C" SEXP _malariavaxTrials_acquire_immunity_cpp(SEXP exposure, SEXP u, SEXP d) {
  BEGIN_CPP11
    return cpp11::as_sexp(acquire_immunity_cpp(cpp11::as_cpp<cpp11::decay_t<doubles>>(exposure), cpp11::as_cpp<cpp11::decay_t<double>>(u), cpp11::as_cpp<cpp11::decay_t<double>>(d)));
  END_CPP11
}
// get_b.cpp
doubles get_b_cpp(doubles ib, double b0, double b1, double ib0, double kb);
extern "C" SEXP _malariavaxTrials_get_b_cpp(SEXP ib, SEXP b0, SEXP b1, SEXP ib0, SEXP kb) {
  BEGIN_CPP11
    return cpp11::as_sexp(get_b_cpp(cpp11::as_cpp<cpp11::decay_t<doubles>>(ib), cpp11::as_cpp<cpp11::decay_t<double>>(b0), cpp11::as_cpp<cpp11::decay_t<double>>(b1), cpp11::as_cpp<cpp11::decay_t<double>>(ib0), cpp11::as_cpp<cpp11::decay_t<double>>(kb)));
  END_CPP11
}
// get_phi.cpp
doubles get_phi_cpp(doubles ica, doubles icm, double phi0, double phi1, double ic0, double kc);
extern "C" SEXP _malariavaxTrials_get_phi_cpp(SEXP ica, SEXP icm, SEXP phi0, SEXP phi1, SEXP ic0, SEXP kc) {
  BEGIN_CPP11
    return cpp11::as_sexp(get_phi_cpp(cpp11::as_cpp<cpp11::decay_t<doubles>>(ica), cpp11::as_cpp<cpp11::decay_t<doubles>>(icm), cpp11::as_cpp<cpp11::decay_t<double>>(phi0), cpp11::as_cpp<cpp11::decay_t<double>>(phi1), cpp11::as_cpp<cpp11::decay_t<double>>(ic0), cpp11::as_cpp<cpp11::decay_t<double>>(kc)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_malariavaxTrials_acquire_immunity_cpp", (DL_FUNC) &_malariavaxTrials_acquire_immunity_cpp, 3},
    {"_malariavaxTrials_get_b_cpp",            (DL_FUNC) &_malariavaxTrials_get_b_cpp,            5},
    {"_malariavaxTrials_get_phi_cpp",          (DL_FUNC) &_malariavaxTrials_get_phi_cpp,          6},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_malariavaxTrials(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
