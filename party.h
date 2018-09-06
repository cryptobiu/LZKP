//
// Created by roee on 9/6/18.
//

#ifndef LZKP_PARTY_H
#define LZKP_PARTY_H


#include <sys/socket.h>
#include <sys/uio.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>


namespace lzkp {


class Party {
public:
  virtual int init(int argc, char **argv) = 0;
  virtual bool runOnline() = 0;
};


}
#endif //LZKP_PARTY_H
