//
// Created by roee on 9/6/18.
//

#ifndef __LZKP_VERIFIER_PARTY_H_FILE__
#define __LZKP_VERIFIER_PARTY_H_FILE__


#include "party.h"

#include <iostream>
#include <cstring> // For memset


namespace lzkp {


template <class FieldType>
class VerifierParty : public Party {
public:
  VerifierParty() : Party() { }
  ~VerifierParty() {
    debug("Closing channel" << std::endl);
    close(sock_);
  }

//  virtual int init(int argc, char **argv);
//  virtual bool runOnline();

protected:
    virtual int initCommunication();

protected:
  std::string ip_;
  int port_;
};

template<class FieldType>
int VerifierParty<FieldType>::initCommunication() {
  struct sockaddr_in serv_addr;

  debug("Initializing communication channel..." << std::endl);
  debug("\tCreating socket... ");

  if ((sock_ = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
    debug("error (" << sock_ << ")" << std::endl);

    return -1;
  }

  std::memset(&serv_addr, '0', sizeof(serv_addr));

  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = inet_addr(ip_.c_str());
  serv_addr.sin_port = htons(port_);

  // Convert IPv4 and IPv6 addresses from text to binary form
  debug("done (" << sock_ << ")" << std::endl);
  debug("\tConverting IP (" << this->ip_.c_str() << ") from text to binary form... ");

  if (int r = inet_pton(AF_INET, this->ip_.c_str(), &serv_addr.sin_addr) <= 0) {
    debug("error (" << r << ")" << std::endl);

    return -2;
  }

  debug("done" << std::endl);
  debug("\tConnecting to " << ip_ << ":" << port_ << "... ");

  if (int r = connect(sock_, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
    debug("error (" << r << ")" << std::endl);

    return -3;
  }

  debug("Initializing communication channel... done" << std::endl);

  return 0;
}


}


#endif // __LZKP_VERIFIER_PARTY_H_FILE__
