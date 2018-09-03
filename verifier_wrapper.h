//
// Created by roee on 9/3/18.
//

#ifndef LZKP_VERIFIER_WRAPPER_H
#define LZKP_VERIFIER_WRAPPER_H

#include <NTL/ZZ_p.h>
#include <NTL/matrix.h>

#include "settings.h"
#include "verifier.h"


namespace lzkp {


class VerifierWrapper {
public:
  VerifierWrapper(const Settings &s, const NTL::Mat<NTL::ZZ_p> &a, const NTL::Vec<NTL::ZZ_p> &t);
  ~VerifierWrapper();

  void r2(const block &h_gamma, std::vector<bool> &E);
  void r4(const std::vector<block> &seed, const std::vector<block> &omegaN, const block &h_pi, block &seed_ell);
  void r6(const block &h_psi, std::vector<int> &i_bar);
  bool r8(const block &seed_e_bar, const std::vector<std::vector<block>> &seed_tree, const std::vector<block> &gamma_i_bar,
          const std::vector<std::vector<NTL::ZZ_p>> &alpha_i_bar, const std::vector<NTL::ZZ_p> &o_i_bar, const std::vector<std::vector<NTL::ZZ_p>> &b_square,
          const std::vector<std::vector<NTL::ZZ_p>> &s);

private:
public:
  // Public known values
  NTL::Mat<NTL::ZZ_p> a_;
  NTL::Vec<NTL::ZZ_p> t_;

  const Settings &set_;
  const int M;
  const int N;
  const uint64_t q;
  const int m;
  const int n;
  const int tau;

  std::vector<Verifier *> verifiers_;

  block h_gamma_;

  std::vector<bool> E_;

//  std::vector<block> seed_;
//  std::vector<block> omegaN_;
  block h_pi_;
  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;

  block h_psi_;

  block seed_e_bar_;
  osuCrypto::PRNG prng_e_bar_;

};


}


#endif //LZKP_VERIFIER_WRAPPER_H
