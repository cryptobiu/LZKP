//
// Created by roee on 9/5/18.
//

#ifndef LZKP_RING31_H
#define LZKP_RING31_H


#include <iostream>


struct Ring31 {
public: //TODO return to private after tesing

  static const unsigned int p = 31;
  unsigned int elem;

public:

  Ring31() { elem = 0; }
//  Ring31(int elem)
//  {
//    this->elem = elem % p;
//  }
  Ring31(uint64_t elem)
  {
    this->elem = elem % p;
  }

  Ring31& operator=(const Ring31 &other) { elem = other.elem; return *this; }
  bool operator==(const Ring31 &other) const { return (other.elem == elem); }
  bool operator!=(const Ring31 &other) const { return !(other.elem == elem); }

  Ring31 operator+(const Ring31 &f2) const
  {
    return Ring31((elem + f2.elem) % p);
  }

  Ring31 operator-(const Ring31 &f2) const
  {
    return Ring31((elem - f2.elem + p) % p);
  }

  Ring31 operator/(const Ring31 &f2) const
  {
    return Ring31((8 * elem) % p);
  }

  Ring31 operator*(const Ring31 &f2) const
  {
    long multLong = ((long)elem * (long)f2.elem) % p;

    return Ring31((int)multLong);
  }

  Ring31& operator+=(const Ring31 &f2) { elem = (f2.elem + elem) % p; return *this; }
  Ring31& operator-=(const Ring31 &f2) { elem = (elem - f2.elem + p) % p; return *this; }
  Ring31& operator*=(const Ring31 &f2)
  {
    long multLong = ((long)elem * (long) f2.elem) % p;

    elem = (int)multLong;

    return *this;
  }
};

inline std::ostream& operator<<(std::ostream& s, const Ring31 &a) { return s << a.elem; }


#endif //LZKP_RING31_H
