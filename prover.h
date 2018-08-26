#ifndef __LZKP_PROVER_H_FILE___
#define __LZKP_PROVER_H_FILE___

#include <NTL/ZZ_p.h>
#include <NTL/matrix.h>

#include "seedtree.h"
#include "settings.h"

namespace lzkp {


class Prover {
public:
  Prover(const Settings &s);
  ~Prover();

  block r1();
  void r3(std::vector<bool> E, std::vector<block> &seed, std::vector<block> &omega, block &h_pi);
  block r5(const std::vector<std::vector<NTL::ZZ_p>> &coefficients);
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
  std::vector<NTL::ZZ_p> secret_;

  std::vector<block> master_seed_;
  std::vector<SeedTree> seed_tree_;
  std::vector<std::vector<block>> r_;
  std::vector<std::vector<std::vector<NTL::ZZ_p>>> b_;
  std::vector<std::vector<std::vector<NTL::ZZ_p>>> b_square_;

  block seed_e_bar_;
  osuCrypto::PRNG prng_e_bar_;

  std::vector<std::vector<block>> gamma_;
  std::vector<block> h_;
  block h_gamma_;

  std::vector<std::vector<std::vector<NTL::ZZ_p>>> s_;
  std::vector<std::vector<std::vector<NTL::ZZ_p>>> alpha_;

  std::vector<bool> E_;
  std::vector<block> g_;
  std::vector<block> gN_; // g_{e,N}
  block h_pi_;

  std::vector<std::vector<NTL::ZZ_p>> o_;
  std::vector<block> w_;
  std::vector<block> psi_;
  block h_psi_;
};


}
#endif