//
// Created by roee on 7/26/18.
//

#ifndef LZKP_VERIFIER_H
#define LZKP_VERIFIER_H

#include <vector>

#include "seedtree.h"
#include "settings.h"

#include <NTL/ZZ_p.h>

namespace lzkp {


class Verifier {
public:
  Verifier(const Settings &s);

  const std::vector<bool> &r2(const block &h_gamma);
  void r4(const std::vector<block> &seed, const std::vector<block> &omega, const block &h_pi, std::vector<std::vector<NTL::ZZ_p>> &coefficients);
  void r6(const block &h_psi, std::vector<int> &i_bar);

private:
public:
  const int M;
  const int N;
  const int q;
  const int m;
  const int n;
  const int tau;

  std::vector<bool> E_;
  block h_gamma_;

  std::vector<block> seed_;
  std::vector<block> omega_;
  block h_pi_;

  std::vector<SeedTree> seed_tree_;
  std::vector<std::vector<block>> r_;
  std::vector<std::vector<std::vector<NTL::ZZ_p>>> b_;
  std::vector<std::vector<std::vector<NTL::ZZ_p>>> b_square_;

  std::vector<std::vector<block>> gamma_;
  std::vector<block> h_;

  block h_psi_;
};


}

#endif //LZKP_VERIFIER_H
