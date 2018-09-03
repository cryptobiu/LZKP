#ifndef LZKP_VERIFIER_H
#define LZKP_VERIFIER_H

#include <NTL/ZZ_p.h>
#include <NTL/matrix.h>

#include "seedtree.h"
#include "settings.h"


namespace lzkp {


class Verifier {
public:
  Verifier(const Settings &s, const NTL::Mat<NTL::ZZ_p> &a, const NTL::Vec<NTL::ZZ_p> &t);
  ~Verifier();

  void r4(); // Local variable seed_ must be set before calling this method
//  bool r8(const block &seed_e_bar, const std::vector<std::vector<block>> &seed_tree, const std::vector<block> &gamma_i_bar,
//          const std::vector<std::vector<NTL::ZZ_p>> &alpha_i_bar, const std::vector<NTL::ZZ_p> &o_i_bar, const std::vector<std::vector<NTL::ZZ_p>> &b_square,
//          const std::vector<std::vector<NTL::ZZ_p>> &s, std::vector<std::vector<block>> &partial_seeds);
  bool r8(const std::vector<block> &seed_tree, const block &gamma_i_bar,
          const std::vector<NTL::ZZ_p> &alpha_i_bar, const NTL::ZZ_p &o_i_bar, const std::vector<NTL::ZZ_p> &b_square,
          const std::vector<NTL::ZZ_p> &s); // Local variables omegaN_ (from round 4), g_, w_,  must be set before calling this method
private:
public:
  const int N;
  const uint64_t q;
  const int m;
  const int n;

  // Public known values
  const NTL::Mat<NTL::ZZ_p> &a_;
  const NTL::Vec<NTL::ZZ_p> &t_;

  // Local values
  block seed_;
  block omegaN_;
  SeedTree seed_tree_;
  std::vector<block> r_;
  std::vector<std::vector<NTL::ZZ_p>> b_;
  std::vector<std::vector<NTL::ZZ_p>> b_square_;
  std::vector<block> gamma_;
  block h_;
  std::vector<NTL::ZZ_p> coefficients_;

  int i_bar_;

  bool reject_;
  std::vector<block> partial_seeds_;
  block pi_;
  block g_;
  std::vector<NTL::ZZ_p> o_;
  NTL::ZZ_p sigma_o;
  block w_;
  block psi_;
};


}
#endif //LZKP_VERIFIER_H
