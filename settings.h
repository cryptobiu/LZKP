//
// Created by roee on 8/15/18.
//

#ifndef LZKP_SETTINGS_H
#define LZKP_SETTINGS_H


class Settings {
public:
  Settings(int M, int N, int q, int m, int t) : M(M), N(N), q(q), m(m), t(t) {}

  const int M;
  const int N;
  const int q;
  const int m;
  const int t;
//  const int M = 5;
//  const int N = 8;
//  const int q = 31;
//  const int m = 10;
};


#endif //LZKP_SETTINGS_H
