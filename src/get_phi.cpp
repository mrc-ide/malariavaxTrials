#include "cpp11.hpp"
using namespace cpp11;

[[cpp11::register]]
doubles get_phi_cpp(doubles ica, doubles icm, double phi0, double phi1, double ic0, double kc) {

  int n = ica.size();

  writable::doubles out(n);

  for(int i = 0; i < n; ++i) {
    out[i] = phi0 * (phi1 + ((1 - phi1) / (1 + pow(((ica[i] + icm[i]) / ic0), kc))));
  }

  return out;

}

