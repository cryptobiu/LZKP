#ifndef __LZKP_PROVER_H_FILE___
#define __LZKP_PROVER_H_FILE___

#include "seedtree.h"
#include "settings.h"

#include <cryptoTools/Crypto/sha1.h>

#include <NTL/ZZ_p.h>

namespace lzkp {


class Prover {
public:
  Prover(const Settings &s);
  ~Prover();

  block r1();
  void r3(std::vector<bool> E);

private:
  const int M;
  const int N;
  const int q;
  const int m;
  const int t;

  std::vector<block> master_seed_;
  std::vector<SeedTree> st_;
  std::vector<std::vector<block>> r_;
  std::vector<std::vector<std::vector<NTL::ZZ>>> b_;
  std::vector<std::vector<std::vector<NTL::ZZ>>> b_square_;

  std::vector<std::vector<block>> gamma_;
};


}
#endif