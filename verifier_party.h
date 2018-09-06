//
// Created by roee on 9/6/18.
//

#ifndef LZKP_VERIFIER_PARTY_H
#define LZKP_VERIFIER_PARTY_H


#include "party.h"

#include <iostream>
#include <cstring>


namespace lzkp {


template <class FieldType>
class VerifierParty : public Party {
public:
  VerifierParty() : sock(0), port(0) { }
  ~VerifierParty() { close(sock); }

  virtual int init(int argc, char **argv);
  virtual bool runOnline();

private:
  int sock;
};

template<class FieldType>
int VerifierParty<FieldType>::init(int argc, char **argv) {
  int port;
  struct sockaddr_in serv_addr;

  if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    return -1;

  std::memset(&serv_addr, '0', sizeof(serv_addr));

  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = inet_addr(ip.c_str());
  serv_addr.sin_port = htons(port);

  // Convert IPv4 and IPv6 addresses from text to binary form
  if(inet_pton(AF_INET, ip.c_str(), &serv_addr.sin_addr)<=0)
    return -2;

  if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
    return -3;

  return 0;
}

template<class FieldType>
bool VerifierParty<FieldType>::runOnline() {}


}


#endif //LZKP_VERIFIER_PARTY_H
