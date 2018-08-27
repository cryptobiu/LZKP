//
// Created by roee on 8/15/18.
//

#ifndef LZKP_SETTINGS_H
#define LZKP_SETTINGS_H

class Settings {
public:
  Settings(int M, int N, uint64_t q, int m, int n, int tau) : M(M), N(N), q(q), m(m), n(n), tau(tau) { }

  const int M;
  const int N;
  const uint64_t q;
  const int m;
  const int n; // n < m
  const int tau;
//  const int M = 5;
//  const int N = 8;
//  const int q = 31;
//  const int m = 10;
};


#endif //LZKP_SETTINGS_H
