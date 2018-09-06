//
// Created by roee on 9/6/18.
//

#ifndef LZKP_PROVER_PARTY_H
#define LZKP_PROVER_PARTY_H


#include "party.h"


namespace lzkp {


template <class FieldType>
class ProverParty : public Party {
public:
  ProverParty() : fd(0), port(0), sock(0) { }
  ~ProverParty() { close(sock); }

  virtual int init(int argc, char **argv);
  virtual bool runOnline();

private:
  int sock;
};

template<class FieldType>
int ProverParty<FieldType>::init(int argc, char **argv) {
  int port, fd;
  struct sockaddr_in address;
  int opt = 1;
  int addrlen = sizeof(address);

  if ((fd = socket(AF_INET, SOCK_STREAM, 0)) == 0)
    return -1;

  // Forcefully attaching socket to the port 8080
  if (setsockopt(fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt)))
    return -2;

  address.sin_family = AF_INET;
  address.sin_addr.s_addr = INADDR_ANY;
  address.sin_port = htons(port);

  // Forcefully attaching socket to the port
  if (bind(fd, (struct sockaddr *)&address, sizeof(address)) < 0)
    return -3;

  if (listen(fd, 3) < 0)
    return -4;

  if ((sock = accept(fd, (struct sockaddr *)&address, (socklen_t*)&addrlen)) < 0)
    return -5;

  return 0;
}

template<class FieldType>
bool ProverParty<FieldType>::runOnline() {

}


}

#endif //LZKP_PROVER_PARTY_H
