//
// Created by roee on 9/6/18.
//

#ifndef __LZKP_FACTORY_H_FILE__
#define __LZKP_FACTORY_H_FILE__


#include "party.h"
#include "cac_prover_party.h"
#include "cac_verifier_party.h"
#include "fields/mersenne.h"
#include "fields/field_15_bit.h"


namespace lzkp {


class Factory {
public:
  Factory() = delete;
  Factory(const Factory&) = delete;
  Factory& operator=(const Factory&) = delete;

  static Party* create(int protocol_type, bool is_prover, int q) {
    if (protocol_type == 0) { // Cut-and-Choose
      if (is_prover) { // Create prover
        if (q == 15) {
          return new CacProverParty<Field15Bit>();
        }
        else if (q == 31) {
          return new CacProverParty<ZpMersenneIntElement>();
        }
//        else if (q == 61) {
//          return new CacProverParty<ZpMersenneLongElement>();
//        }
      }
      else { // Create verifier
        if (q == 15) {
          return new CacVerifierParty<Field15Bit>();
        }
        else if (q == 31) {
          return new CacVerifierParty<ZpMersenneIntElement>();
        }
//        else if (q == 61) {
//          return new CacVerifierParty<ZpMersenneLongElement>();
//        }
      }
    }

    return nullptr;
  }
};


}


#endif // __LZKP_FACTORY_H_FILE__
