#ifndef LZKP_SETTINGS_H
#define LZKP_SETTINGS_H


namespace lzkp {


class Settings {
public:
  Settings(int M, int N, uint64_t q, int m, int n, int tau) : M(M), N(N), q(q), m(m), n(n), tau(tau) {}

  const int M;
  const int N;
  const uint64_t q;
  const int m;
  const int n; // n < m
  const int tau;
};


}
#endif //LZKP_SETTINGS_H
