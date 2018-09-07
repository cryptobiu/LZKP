//
// Created by roee on 9/6/18.
//

#ifndef LZKP_PROVER_PARTY_H
#define LZKP_PROVER_PARTY_H


#include "party.h"

#include <iostream>

#include <cryptoTools/Crypto/PRNG.h>


namespace lzkp {


template <class FieldType>
class ProverParty : public Party {
public:
  ProverParty() : Party() { }
  ~ProverParty() {
#ifdef DEBUG
std::cout << "Closing channel" << std::endl;
#endif
  close(sock_);
  }

protected:
  virtual int initCommunication();
  virtual int generateData() = 0;

protected:
  int port_;
};

template<class FieldType>
int ProverParty<FieldType>::initCommunication() {
  int fd;
  struct sockaddr_in address;
  int opt = 1;
  int addrlen = sizeof(address);

#ifdef DEBUG
  std::cout << "Initializing communication channel..." << std::endl;
  std::cout << "\tCreating socket... ";
#endif

  if ((fd = socket(AF_INET, SOCK_STREAM, 0)) == 0) {
#ifdef DEBUG
    std::cout << "error (" << fd << ")" << std::endl;
#endif
    return -1;
  }

#ifdef DEBUG
  std::cout << "done (" << fd << ")" << std::endl;
  std::cout << "\tSetting socket options... ";
#endif

  // Forcefully attaching socket to the port 8080
  if (int r = setsockopt(fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt))) {
#ifdef DEBUG
    std::cout << "error (" << r << ")" << std::endl;
#endif
    return -2;
  }

  address.sin_family = AF_INET;
  address.sin_addr.s_addr = INADDR_ANY;
  address.sin_port = htons(port_);

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\tBinding socket... ";
#endif

  // Forcefully attaching socket to the port
  if (int r = bind(fd, (struct sockaddr *)&address, sizeof(address)) < 0) {
#ifdef DEBUG
    std::cout << "error (" << r << ")" << std::endl;
#endif
    return -3;
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\tListening on port " << port_ << "... ";
#endif
  if (int r = listen(fd, 3) < 0) {
#ifdef DEBUG
    std::cout << "error (" << r << ")" << std::endl;
#endif
    return -4;
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\tAccepting connction... ";
#endif
  if ((sock_ = accept(fd, (struct sockaddr *)&address, (socklen_t*)&addrlen)) < 0) {
#ifdef DEBUG
    std::cout << "error (" << sock_ << ")" << std::endl;
#endif
    return -5;
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "Initializing communication channel... done" << std::endl;
#endif

  return 0;
}


}

#endif //LZKP_PROVER_PARTY_H
