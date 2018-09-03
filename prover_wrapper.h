//
// Created by roee on 9/3/18.
//

#ifndef LZKP_PROVER_WRAPPER_H
#define LZKP_PROVER_WRAPPER_H

#include <NTL/ZZ_p.h>
#include <NTL/matrix.h>

#include "settings.h"
#include "prover.h"

namespace lzkp {


class ProverWrapper {
public:
  ProverWrapper(const Settings &s, const NTL::Mat<NTL::ZZ_p> &a, const NTL::Vec<NTL::ZZ_p> &t, const NTL::Vec<NTL::ZZ_p> &secret);
  ~ProverWrapper();

  void r1(block &h_gamma);
  void r3(const std::vector<bool> &E, std::vector<block> &seed, std::vector<block> &omegaN, block &h_pi);
  void r5(const block &seed_ell, block &h_psi);
  void r7(const std::vector<int> &i_bar, block &seed_e_bar, std::vector<std::vector<block>> &seed_tree,
          std::vector<block> &gamma_i_bar, std::vector<std::vector<NTL::ZZ_p>> &alpha_i_bar, std::vector<NTL::ZZ_p> &o_i_bar,
          std::vector<std::vector<NTL::ZZ_p>> &b_square, std::vector<std::vector<NTL::ZZ_p>> &s);

private:
public:
  // Public known values
  NTL::Mat<NTL::ZZ_p> a_;
  NTL::Vec<NTL::ZZ_p> t_;

  // Prover's secret
  NTL::Vec<NTL::ZZ_p> secret_;

  const Settings &set_;
  const int M;
  const int N;
  const uint64_t q;
  const int m;
  const int n;
  const int tau;

  std::vector<Prover *> provers_;

  std::vector<block> master_seed_;
  block h_gamma_;

  std::vector<bool> E_;
  block seed_e_bar_;
  osuCrypto::PRNG prng_e_bar_;
  std::vector<block> g_;
  block h_pi_;

  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;
  block h_psi_;
};


}


#endif //LZKP_PROVER_WRAPPER_H
