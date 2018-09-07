//
// Created by roee on 9/6/18.
//

#ifndef LZKP_VERIFIER_PARTY_H
#define LZKP_VERIFIER_PARTY_H


#include "party.h"

#include <iostream>
#include <cstring> // For memset


namespace lzkp {


template <class FieldType>
class VerifierParty : public Party {
public:
  VerifierParty() : Party() { }
  ~VerifierParty() {
#ifdef DEBUG
    std::cout << "Closing channel" << std::endl;
#endif
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

#ifdef DEBUG
  std::cout << "Initializing communication channel..." << std::endl;
  std::cout << "\tCreating socket... ";
#endif

  if ((sock_ = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
#ifdef DEBUG
  std::cout << "error (" << sock_ << ")" << std::endl;
#endif
    return -1;
  }

  std::memset(&serv_addr, '0', sizeof(serv_addr));

  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = inet_addr(ip_.c_str());
  serv_addr.sin_port = htons(port_);

  // Convert IPv4 and IPv6 addresses from text to binary form
#ifdef DEBUG
  std::cout << "done (" << sock_ << ")" << std::endl;
  std::cout << "\tconverting IP (" << this->ip_.c_str() << ") from text to binary form... ";
#endif

  if (int r = inet_pton(AF_INET, this->ip_.c_str(), &serv_addr.sin_addr) <= 0) {
#ifdef DEBUG
    std::cout << "error (" << r << ")" << std::endl;
#endif
    return -2;
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\tconnecting to " << ip_ << ":" << port_ << "... ";
#endif

  if (int r = connect(sock_, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
#ifdef DEBUG
    std::cout << "error (" << r << ")" << std::endl;
#endif
    return -3;
  }

#ifdef DEBUG
  std::cout << "Initializing communication channel... done" << std::endl;
#endif

  return 0;
}


}


#endif //LZKP_VERIFIER_PARTY_H
