//
// Created by roee on 9/6/18.
//

#ifndef __LZKP_PROVER_PARTY_H_FILE__
#define __LZKP_PROVER_PARTY_H_FILE__


#include "party.h"

#include <iostream>

#include <cryptoTools/Crypto/PRNG.h>


namespace lzkp {


template <class FieldType>
class ProverParty : public Party {
public:
  ProverParty() : Party() { }
  ~ProverParty() {
    debug("Closing channel" << std::endl);
    close(sock_);
  }

protected:
  int initCommunication();
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

  debug("Initializing communication channel..." << std::endl);
  debug("\tCreating socket... ");

  if ((fd = socket(AF_INET, SOCK_STREAM, 0)) == 0) {
    debug("error (" << fd << ")" << std::endl);

    return -1;
  }

  debug("done (" << fd << ")" << std::endl);
  debug("\tSetting socket options... ");

  // Forcefully attaching socket to the port 8080
  if (int r = setsockopt(fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt))) {
    debug("error (" << r << ")" << std::endl);

    return -2;
  }

  address.sin_family = AF_INET;
  address.sin_addr.s_addr = INADDR_ANY;
  address.sin_port = htons(port_);

  debug("done" << std::endl);
  debug("\tBinding socket... ");

  // Forcefully attaching socket to the port
  if (int r = bind(fd, (struct sockaddr *)&address, sizeof(address)) < 0) {
    debug("error (" << r << ")" << std::endl);

    return -3;
  }

  debug("done" << std::endl);
  debug("\tListening on port " << port_ << "... ");

  if (int r = listen(fd, 3) < 0) {
    debug("error (" << r << ")" << std::endl);

    return -4;
  }

  debug("done" << std::endl);
  debug("\tAccepting connction... ");

  if ((sock_ = accept(fd, (struct sockaddr *)&address, (socklen_t*)&addrlen)) < 0) {
    debug("error (" << sock_ << ")" << std::endl);

    return -5;
  }

  debug("done" << std::endl);
  debug("Initializing communication channel... done" << std::endl);

  return 0;
}


}

#endif // __LZKP_PROVER_PARTY_H_FILE__
