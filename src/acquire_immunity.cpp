#include "cpp11.hpp"
using namespace cpp11;

[[cpp11::register]]
doubles acquire_immunity_cpp(doubles exposure, double u, double d) {
  int n = exposure.size();

  writable::doubles immunity(n);
  immunity[0] = 0;

  for(int i = 1; i < n; ++i) {
    immunity[i] = immunity[i - 1] + (exposure[i] / (exposure[i] * u + 1)) - (immunity[i - 1] / d);
  }

  return immunity;
}
