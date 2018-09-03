#ifndef __LZKP_PROVER_H_FILE___
#define __LZKP_PROVER_H_FILE___

#include <NTL/ZZ_p.h>
#include <NTL/matrix.h>

#include "seedtree.h"
#include "settings.h"


namespace lzkp {


class Prover {
public:
  Prover(const Settings &s, const NTL::Mat<NTL::ZZ_p> &a, const NTL::Vec<NTL::ZZ_p> &t, const NTL::Vec<NTL::ZZ_p> &secret);
  ~Prover();

  void r1(block &h_gamma);
  void r3(const std::vector<bool> &E, std::vector<block> &seed, std::vector<block> &omegaN, block &h_pi);
  void r5(const block &seed_ell, block &h_psi);
  void r7(const std::vector<int> &i_bar, block &seed_e_bar, std::vector<std::vector<block>> &seed_tree,
          std::vector<block> &gamma_i_bar, std::vector<std::vector<NTL::ZZ_p>> &alpha_i_bar, std::vector<NTL::ZZ_p> &o_i_bar,
          std::vector<std::vector<NTL::ZZ_p>> &b_square, std::vector<std::vector<NTL::ZZ_p>> &s);

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

  // Prover's secret
  NTL::Vec<NTL::ZZ_p> secret_;

  // Local values
  std::vector<block> master_seed_;
  std::vector<SeedTree> seed_tree_;
  std::vector<std::vector<block>> r_;
  std::vector<std::vector<std::vector<NTL::ZZ_p>>> b_;
  std::vector<std::vector<std::vector<NTL::ZZ_p>>> b_square_;
  std::vector<std::vector<block>> gamma_;
  std::vector<block> h_;
  block h_gamma_;

  std::vector<bool> E_;

  block seed_e_bar_;
  osuCrypto::PRNG prng_e_bar_;
  std::vector<block> g_;
  std::vector<std::vector<std::vector<NTL::ZZ_p>>> s_;
  std::vector<block> gN_; // g_{e,N}
  std::vector<std::vector<std::vector<NTL::ZZ_p>>> alpha_;
  std::vector<std::vector<NTL::ZZ_p>> alpha_sum_;
  std::vector<block> seed_;
  std::vector<block> omegaN_;
  block h_pi_;

  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;
  std::vector<std::vector<NTL::ZZ_p>> coefficients_;
  std::vector<std::vector<NTL::ZZ_p>> o_;
  std::vector<block> w_;
  std::vector<block> psi_;
  block h_psi_;
};


}
#endif