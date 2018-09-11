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
#include <sys/ioctl.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#ifdef DEBUG
#include <boost/type_index.hpp>
#endif

#include "block.h"


namespace lzkp {


class Party {
public:
  Party() : sock_(0) { }
  virtual ~Party() { };

  virtual int init(int argc, const char* const argv[]) = 0;
  virtual bool runOnline() = 0;

protected:
  virtual int parseArguments(int argc, const char *const *argv) = 0;
  virtual int initCommunication() = 0;
  virtual int negotiateParameters() = 0;

  int sock_;
};


}
#endif //LZKP_PARTY_H
