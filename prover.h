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

  void r1(const block &master_seed);
  void r3(); // Local variable g_ must be set before calling this method
  void r5(); // Local variable coefficients_ must be set before calling this method
  void r7(std::vector<block> &seed_tree, block &gamma_i_bar, std::vector<NTL::ZZ_p> &alpha_i_bar,
          NTL::ZZ_p &o_i_bar, std::vector<NTL::ZZ_p> &b_square, std::vector<NTL::ZZ_p> &s);
             // Local variable i_bar_ must be set before calling this method

private:
public:
  const int N;
  const uint64_t q;
  const int m;
  const int n;

  // Public known values
  const NTL::Mat<NTL::ZZ_p> &a_;
  const NTL::Vec<NTL::ZZ_p> &t_;

  // Prover's secret
  const NTL::Vec<NTL::ZZ_p> &secret_;

  // Local values
  block master_seed_;
  SeedTree seed_tree_;
  std::vector<block> r_;
  std::vector<std::vector<NTL::ZZ_p>> b_;
  std::vector<std::vector<NTL::ZZ_p>> b_square_;
  std::vector<block> gamma_;
  block h_;

  block g_;
  std::vector<std::vector<NTL::ZZ_p>> s_;
  std::vector<std::vector<NTL::ZZ_p>> alpha_;
  std::vector<NTL::ZZ_p> alpha_sum_;
  block gN_; // g_{e,N}
  block omegaN_;
  block pi_;

  std::vector<NTL::ZZ_p> coefficients_;
  std::vector<NTL::ZZ_p> o_;
  block w_;
  block psi_;

  int i_bar_;
};


}
#endif