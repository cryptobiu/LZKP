#ifndef LZKP_SETTINGS_H
#define LZKP_SETTINGS_H

#include <cstdint>


namespace lzkp {


class Settings {
public:
  Settings(int M, int tau, int N, uint64_t q, int n, int m) : M(M), tau(tau), N(N), n(n), m(m) {}

  const int M;
  const int tau;
  const int N;
//  const uint64_t q;
  const int n; // n < m
  const int m;
};


}
#endif //LZKP_SETTINGS_H
