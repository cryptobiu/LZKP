//
// Created by roee on 9/6/18.
//

#ifndef __LZKP_PARTY_H_FILE__
#define __LZKP_PARTY_H_FILE__


#include <sys/socket.h>
#include <sys/uio.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <cerrno>
#include <chrono>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#ifdef DEBUG
  #define debug(x) std::cerr << x
  #include <boost/type_index.hpp>
#else
  #define debug(x)
#endif

#include "block.h"
#include "utils.h"


namespace lzkp {


class Party {
public:
  Party() : sock_(0), RTT_(0), tot_computation_time_(0) { }
  virtual ~Party() { };

  virtual int init(int argc, const char* const argv[]) = 0;
  virtual bool runOnline() = 0;

protected:
  virtual int parseArguments(int argc, const char *const *argv) = 0;
  virtual int initCommunication() = 0;
  virtual int sync() = 0;
  virtual int negotiateParameters() = 0;
  int writeWrapper(void *buf, size_t count);
  int readWrapper(void *buf, size_t count);
  int writevWrapper(iovec *iov, int iovcnt, ssize_t nexpected);
  int readvWrapper(iovec *iov, int iovcnt, ssize_t nexpected);
  void startComputationClock();
  void stopComputationClock();

  int sock_;
  std::chrono::time_point<std::chrono::high_resolution_clock> computation_clock_;

public:
  int RTT_;
  int tot_computation_time_;
  int time_cut_and_choose_;
  int time_eq_1;
};


}
#endif // __LZKP_PARTY_H_FILE__
