//
// Created by roee on 7/26/18.
//

#include "verifier.h"

#include <cryptoTools/Crypto/PRNG.h>

#include "seedtree.h"

using namespace lzkp;

Verifier::Verifier(const Settings &s) : M(s.M), N(s.N), q(s.q), m(s.m), t(s.t) {
  E_.resize(M);
}

const std::vector<bool> &Verifier::r2() {
  osuCrypto::PRNG prng;

  prng.SetSeed(osuCrypto::sysRandomSeed());

  for (auto i = M - t + 1; i <= M; ++i) {
    long r = prng.get<osuCrypto::u64>() % i;

    if (!E_[r])
      E_[r] = true;
    else
      E_[i - 1] = true;
  }

  return E_;
}
