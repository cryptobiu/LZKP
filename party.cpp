//
// Created by lzkp on 9/12/18.
//

#include "party.h"


using namespace lzkp;


int Party::writeWrapper(void *buf, size_t count) {
  ssize_t nwritten;

  nwritten = write(this->sock_, buf, count);

  std::cout << nwritten << std::endl;

  if (nwritten == -1) {
    std::cerr << "writev error (#" << errno << ")" << std::endl;

    return -1;
  }

  assert (nwritten == (ssize_t)count);

  return 0;
}

int Party::readWrapper(void *buf, size_t count) {
  ssize_t nread;
  int nready;

  nready = -1;

  while (nready < (int)count) {
    std::cout << nready << " " << count << std::endl;
    std::cout << ioctl(this->sock_, FIONREAD, &nready) << std::endl;
  }

  nread = read(this->sock_, buf, count);

  if (nread == -1) {
    std::cerr << "readv error (#" << errno << ")" << std::endl;

    return -1;
  }

  assert (nread == (ssize_t)count);

  return 0;
}

int Party::writevWrapper(const iovec *iov, int iovcnt, ssize_t nexpected) {
  ssize_t nwritten;

  nwritten = writev(this->sock_, iov, iovcnt);

  if (nwritten == -1) {
    std::cerr << "writev error (#" << errno << ")" << std::endl;

    return -1;
  }

  assert (nwritten == nexpected);

  return 0;
}

int Party::readvWrapper(const iovec *iov, int iovcnt, ssize_t nexpected) {
  ssize_t nread;
  int nready;

  nready = -1;

  while (nready < nexpected) {
    std::cout << nready << " " << nexpected << std::endl;
    ioctl(this->sock_, FIONREAD, &nready);
  }

  nread = readv(this->sock_, iov, iovcnt);

  if (nread == -1) {
    std::cerr << "readv error (#" << errno << ")" << std::endl;

    return -1;
  }

  assert (nread == nexpected);

  return 0;
}