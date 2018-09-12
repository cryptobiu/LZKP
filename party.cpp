//
// Created by lzkp on 9/12/18.
//

#include "party.h"


using namespace lzkp;


int Party::writeWrapper(void *buf, size_t count) {
  ssize_t nwritten;

  nwritten = write(this->sock_, buf, count);

  if (nwritten == -1) {
    std::cerr << "writev error (#" << errno << ")" << std::endl;

    exit(1);
  }

  assert (nwritten == (ssize_t)count);

  return 0;
}

int Party::readWrapper(void *buf, size_t count) {
  ssize_t nread;
  size_t ntot = 0;

  while (count) {
    nread = read(this->sock_, (char *)buf + ntot, count);

    if (nread == -1) {
      std::cerr << "readv error (#" << errno << ")" << std::endl;

      exit(1);
    }

    ntot += nread;
    count -= nread;
  }

  return 0;
}

int Party::writevWrapper(iovec *iov, int iovcnt, ssize_t nexpected) {
  ssize_t nwritten;

  nwritten = writev(this->sock_, iov, iovcnt);

  if (nwritten == -1) {
    std::cerr << "writev error (#" << errno << ")" << std::endl;

    exit(1);
  }

  assert (nwritten == nexpected);

  return 0;
}

int Party::readvWrapper(iovec *iov, int iovcnt, ssize_t nexpected) {
  ssize_t nread;

  while (nexpected) {
    nread = readv(this->sock_, iov, iovcnt);

    if (nread == -1) {
      std::cerr << "readv error (#" << errno << ")" << std::endl;

      exit(1);
    }

    nexpected -= nread;

    while (nread) {
      if (nread >= (int)iov->iov_len) {
        nread -= iov->iov_len;
        iov++;
        iovcnt--;

        continue;
      }

      iov->iov_base = (char *)iov->iov_base + nread;
      iov->iov_len -= nread;

      break;
    }
  }

  return 0;
}