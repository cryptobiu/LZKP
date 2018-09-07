#ifndef __LZKP_SETTINGS_H_FILE__
#define __LZKP_SETTINGS_H_FILE__


namespace lzkp {


struct Settings {
public:
  Settings() : M(0), tau(0), N(0), n(0), m(0) { }
  Settings(int M, int tau, int N, int n, int m) : M(M), tau(tau), N(N), n(n), m(m) { }

  int M; // M > 0
  int tau; // tau <= M
  int N;
  int n; // n < m
  int m;

//  const Settings& operator=(const Settings &rhs) {
//    M = rhs.M;
//    tau = rhs.tau;
//    N = rhs.N;
//    n = rhs.n;
//    m = rhs.m;
//
//    return *this;
//  }
};


}


#endif // __LZKP_SETTINGS_H_FILE__
