//
// Created by lzkp on 9/7/18.
//

#ifndef __LZKP_FIELD_15_BIT_H_FILE__
#define __LZKP_FIELD_15_BIT_H_FILE__


#include <iostream>


struct Field15Bit {
//private:
public: //TODO return to private after tesing
  static const uint16_t p = 0x7FED; // 2^15 - 19
  uint16_t elem;

public:
  Field15Bit() { elem = 0; }
  Field15Bit(uint64_t elem) {
    this->elem = elem % p;
  }

  Field15Bit& operator=(const Field15Bit &other) { elem = other.elem; return *this; }
  bool operator==(const Field15Bit &other) const { return (other.elem == elem); }
  bool operator!=(const Field15Bit &other) const { return !(other.elem == elem); }

  Field15Bit operator+(const Field15Bit &f2) const {
    uint32_t temp = (uint32_t)elem + (uint32_t)f2.elem;

    if (temp >= p) {
      temp -= p;
    }

    return Field15Bit(temp);
  }

  Field15Bit operator-(const Field15Bit &f2) const {
    int32_t temp = (int32_t)elem - (int32_t)f2.elem;

    if (temp < 0) {
      return Field15Bit(temp + p);
    }
    else {
      return Field15Bit(temp);
    }
  }

  Field15Bit operator/(const Field15Bit &f2) const {
    // Code taken from NTL for the function XGCD
    int a = f2.elem;
    int b = p;
    long s;

    int  u, v, q, r;
    long u0, v0, u1, v1, u2, v2;

    int aneg = 0;

    if (a < 0) {
      a = -a;
      aneg = 1;
    }

    if (b < 0) {
      b = -b;
    }

    u1=1; v1=0;
    u2=0; v2=1;
    u = a; v = b;

    while (v != 0) {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 =  u1 - q*u2;
      v2 = v1- q*v2;
      u1 = u0;
      v1 = v0;
    }

    if (aneg)
      u1 = -u1;

    s = u1;

    if (s < 0)
      s =  s + p;

    Field15Bit inverse(s);

    return inverse * (*this);
  }

  Field15Bit operator*(const Field15Bit &f2) const {
    return Field15Bit(((uint32_t)elem * (uint32_t)f2.elem) % p);
  }

  Field15Bit& operator+=(const Field15Bit &f2) {
    uint32_t temp = (uint32_t)elem + (uint32_t)f2.elem;

    if (temp >= p) {
      temp -= p;
    }

    elem = temp;

    return *this;
  }

  Field15Bit& operator-=(const Field15Bit &f2) {
    int32_t temp = (int32_t)elem - (int32_t)f2.elem;

    if (temp < 0) {
      elem = temp + p;
    }
    else {
      elem = temp;
    }

    return *this;
  }

  Field15Bit& operator*=(const Field15Bit &f2)
  {
    elem = ((uint32_t)elem * (uint32_t)f2.elem) % p;

    return *this;
  }
};


inline std::ostream& operator<<(std::ostream& s, const Field15Bit &a) { return s << a.elem; }


#endif // __LZKP_FIELD_15_BIT_H_FILE__
