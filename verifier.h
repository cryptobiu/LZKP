//
// Created by roee on 7/26/18.
//

#ifndef LZKP_VERIFIER_H
#define LZKP_VERIFIER_H

#include <vector>

#include "seedtree.h"
#include "settings.h"

#include <NTL/ZZ_p.h>
#include <NTL/matrix.h>

namespace lzkp {


class Verifier {
public:
  Verifier(const Settings &s, const NTL::Mat<NTL::ZZ_p> &a, const NTL::Vec<NTL::ZZ_p> &t);

  void r2(const block &h_gamma, std::vector<bool> &E);
  void r4(const std::vector<block> &seed, const std::vector<block> &omegaN, const block &h_pi, block &seed_ell);
  void r6(const block &h_psi, std::vector<int> &i_bar);
  bool r8(const block &seed_e_bar, const std::vector<std::vector<block>> &seed_tree, const std::vector<block> &gamma_i_bar,
          const std::vector<std::vector<NTL::ZZ_p>> &alpha_i_bar, const std::vector<NTL::ZZ_p> &o_i_bar, const std::vector<std::vector<NTL::ZZ_p>> &b_square,
          const std::vector<std::vector<NTL::ZZ_p>> &s, std::vector<std::vector<block>> &partial_seeds);

private:
public:
  const int M;
  const int N;
  const uint64_t q;
  const int m;
  const int n;
  const int tau;

  // Public known values
  NTL::Mat<NTL::ZZ_p> a_;
  NTL::Vec<NTL::ZZ_p> t_;

  // Local values
  block h_gamma_;
  std::vector<bool> E_;

  std::vector<block> seed_;
  std::vector<block> omegaN_;
  block h_pi_;
  std::vector<SeedTree> seed_tree_;
  std::vector<std::vector<block>> r_;
  std::vector<std::vector<std::vector<NTL::ZZ_p>>> b_;
  std::vector<std::vector<std::vector<NTL::ZZ_p>>> b_square_;
  std::vector<std::vector<block>> gamma_;
  std::vector<block> h_;
  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;
  std::vector<std::vector<NTL::ZZ_p>> coefficients_;

  block h_psi_;
  std::vector<int> i_bar_;

  std::vector<std::vector<block>> partial_seeds_;
  block seed_e_bar_;
  osuCrypto::PRNG prng_e_bar_;
  std::vector<block> g_;
  std::vector<std::vector<NTL::ZZ_p>> o_;
  std::vector<block> w_;
  std::vector<block> psi_;
};


}

#endif //LZKP_VERIFIER_H
