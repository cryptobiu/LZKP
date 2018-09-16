//
// Created by roee on 9/6/18.
//

#ifndef __LZKP_PROVER_PARTY_H_FILE__
#define __LZKP_PROVER_PARTY_H_FILE__


#include "party.h"

#include <iostream>

#include <cryptoTools/Crypto/PRNG.h>


namespace lzkp {


class ProverParty : public Party {
public:
  ProverParty() : Party(), port_(0) { }
  ~ProverParty() {
    debug("Closing channel... ");
    close(sock_);
    debug("done" << std::endl);
  }

protected:
  int initCommunication();
  virtual int sync();
  virtual int generateData() = 0;

protected:
  int port_;
};


}

#endif // __LZKP_PROVER_PARTY_H_FILE__
