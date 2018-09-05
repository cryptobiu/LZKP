//
// Created by roee on 9/5/18.
//

#ifndef LZKP_RINGNTL_H
#define LZKP_RINGNTL_H


#include <NTL/ZZ_p.h>


struct RingNTL {
public: //TODO return to private after tesing

  static const unsigned int p = 31;
  NTL::ZZ_p elem;

public:

  RingNTL() { elem = NTL::ZZ_p(0); }
//  Ring31(int elem)
//  {
//    this->elem = elem % p;
//  }
  RingNTL(uint64_t elem)
  {
    this->elem = NTL::ZZ_p(elem);
  }

  RingNTL(NTL::ZZ_p elem)
  {
    this->elem = elem;
  }

  RingNTL& operator=(const RingNTL &other) { elem = other.elem; return *this; }
  bool operator==(const RingNTL &other) const { return (other.elem == elem); }
  bool operator!=(const RingNTL &other) const { return !(other.elem == elem); }

  RingNTL operator+(const RingNTL &f2) const
  {
    return RingNTL(elem + f2.elem);
  }

  RingNTL operator-(const RingNTL &f2) const
  {
    return RingNTL(elem - f2.elem);
  }

  RingNTL operator/(const RingNTL &f2) const
  {
    return RingNTL(elem / f2.elem);
  }

  RingNTL operator*(const RingNTL &f2) const
  {
    return RingNTL(elem * f2.elem);
  }

  RingNTL& operator+=(const RingNTL &f2) { elem = (f2.elem + elem); return *this; }
  RingNTL& operator-=(const RingNTL &f2) { elem = (elem - f2.elem); return *this; }
  RingNTL& operator*=(const RingNTL &f2)
  {
    this->elem *= f2.elem;

    return *this;
  }
};

inline std::ostream& operator<<(std::ostream& s, const RingNTL &a) { return s << a.elem; }


#endif //LZKP_RINGNTL_H
