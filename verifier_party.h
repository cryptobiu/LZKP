//
// Created by roee on 9/6/18.
//

#ifndef __LZKP_VERIFIER_PARTY_H_FILE__
#define __LZKP_VERIFIER_PARTY_H_FILE__


#include "party.h"

#include <iostream>
#include <cstring> // For memset


namespace lzkp {


class VerifierParty : public Party {
public:
  VerifierParty() : Party(), port_(0), is_accepted_(false) { }
  ~VerifierParty() {
    debug("Closing channel... ");
    close(sock_);
    debug("done" << std::endl);
  }

  bool isAccepted() const { return is_accepted_; }

protected:
  int initCommunication();
  virtual int sync();

protected:
  std::string ip_;
  int port_;
  bool is_accepted_;
};


}


#endif // __LZKP_VERIFIER_PARTY_H_FILE__
