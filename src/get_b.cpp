#include "cpp11.hpp"
using namespace cpp11;

[[cpp11::register]]
doubles get_b_cpp(doubles ib, double b0, double b1, double ib0, double kb) {

  int n = ib.size();

  writable::doubles out(n);

  for(int i = 0; i < n; ++i) {
    out[i] = b0 * (b1 + ((1 - b1) / (1 + pow((ib[i] / ib0), kb))));
  }

  return out;

}
