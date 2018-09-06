//
// Created by roee on 9/6/18.
//

#ifndef LZKP_FACTORY_H
#define LZKP_FACTORY_H

#include "party.h"
#include "prover_party.h"
#include "verifier_party.h"
#include "Mersenne.h"

namespace lzkp {


static class Factory {
public:
  static Party* create(int party_type, int q) {
    if (party_type == 0) { // Create prover
      if (q == 31) {
        return new ProverParty<ZpMersenneIntElement>();
      }
      else if (q == 61) {
        return new ProverParty<ZpMersenneLongElement>();
      }
    }
    else { // Create verifier
      if (q == 31) {
        return new VerifierParty<ZpMersenneIntElement>();
      }
      else if (q == 61) {
        return new VerifierParty<ZpMersenneLongElement>();
      }
    }

    return nullptr;
  }
};


}


#endif //LZKP_FACTORY_H
