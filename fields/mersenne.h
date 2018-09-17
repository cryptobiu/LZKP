//
// Created by moriya on 01/10/17.
// Adjusted by Roee on 04/09/18.
//

#ifndef __LZKP_MERSENNE_H_FILE__
#define __LZKP_MERSENNE_H_FILE__


#include <iostream>
#include <vector>
#include <x86intrin.h>
#include <gmp.h>


typedef uint8_t byte;
typedef unsigned __int128 uint128_t;

struct ZpMersenneIntElement {

//private:
public: //TODO return to private after tesing

  static const unsigned int p = 0x7FFFFFFF; // 2^31 - 1;
  unsigned int elem;

public:

  ZpMersenneIntElement() { elem = 0; }
  ZpMersenneIntElement(int elem)
  {
    this->elem = elem;
    if (this->elem < p) {
      return;
    }
    this->elem -= p;
    if (this->elem < p) {
      return;
    }
    this->elem -= p;
  }
  ZpMersenneIntElement(const ZpMersenneIntElement &rhs) {
    this->elem = rhs.elem;
  }

  ZpMersenneIntElement& operator=(const ZpMersenneIntElement &other) { elem = other.elem; return *this; }
  bool operator==(const ZpMersenneIntElement &other) const { return (other.elem == elem); }
  bool operator!=(const ZpMersenneIntElement &other) const { return !(other.elem == elem); }

  ZpMersenneIntElement operator+(const ZpMersenneIntElement &f2) const
  {
    ZpMersenneIntElement answer;

    answer.elem = (elem + f2.elem);

    if (answer.elem >= p)
      answer.elem -= p;

    return answer;
  }

  ZpMersenneIntElement operator-(const ZpMersenneIntElement &f2) const
  {
    ZpMersenneIntElement answer;

    int temp =  (int)elem - (int)f2.elem;

    if (temp < 0) {
      answer.elem = temp + p;
    }
    else {
      answer.elem = temp;
    }

    return answer;
  }

  ZpMersenneIntElement operator/(const ZpMersenneIntElement &f2) const
  {
    // Code taken from NTL for the function XGCD
    int a = f2.elem;
    int b = p;
    long s;

    int  u, v, q, r;
    long u0, v0, u1, v1, u2, v2;

    int aneg = 0; //, bneg = 0;

    if (a < 0) {
      //if (a < -NTL_MAX_LONG)
      //  Error("XGCD: integer overflow");
      a = -a;
      aneg = 1;
    }

    if (b < 0) {
      //if (b < -NTL_MAX_LONG)
      //  Error("XGCD: integer overflow");
      b = -b;
//      bneg = 1;
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

    ZpMersenneIntElement inverse(s);

    return inverse * (*this);
  }

  ZpMersenneIntElement operator*(const ZpMersenneIntElement &f2) const
  {
    ZpMersenneIntElement answer;

    long multLong = (long)elem * (long)f2.elem;

    // Get the bottom 31 bit
    unsigned int bottom = multLong & p;

    // Get the top 31 bits
    unsigned int top = (multLong >> 31);

    answer.elem = bottom + top;

    // Maximum value of 2p-2
    if (answer.elem >= p)
      answer.elem -= p;

    //return ZpMersenneIntElement((bottom + top) %p);
    return answer;
  }

  ZpMersenneIntElement& operator+=(const ZpMersenneIntElement &f2) { elem = (f2.elem + elem) % p; return *this; }
  ZpMersenneIntElement& operator-=(const ZpMersenneIntElement &f2) { elem = (elem - f2.elem + p) % p; return *this; }
  ZpMersenneIntElement& operator*=(const ZpMersenneIntElement &f2)
  {
    long multLong = (long)elem * (long) f2.elem;

    // Get the bottom 31 bit
    unsigned int bottom = multLong & p;

    // Get the top 31 bits
    unsigned int top = (multLong >> 31) ;

    elem = bottom + top;

    // Maximum value of 2p-2
    if (elem >= p)
      elem -= p;

    return *this;
  }

  static ZpMersenneIntElement dotProdct(ZpMersenneIntElement **&mat_a, const std::vector<std::vector<ZpMersenneIntElement>> &mat_b, int row, int col, int length) {
    uint128_t v = 0;

    for (int k = 0; k < length; ++k) { // OK up to length = 2^68
      v += (uint64_t)mat_a[row][k].elem * (uint64_t)mat_b[k][col].elem;
    }

    return ZpMersenneIntElement(v % p);
  }

  static ZpMersenneIntElement dotProdct(ZpMersenneIntElement **&mat_a, const std::vector<ZpMersenneIntElement> &vec_b, int row, int col, int length) {
    uint128_t v = 0;

    for (int k = 0; k < length; ++k) { // OK up to length = 2^68
      v += (uint64_t)mat_a[row][k].elem * (uint64_t)vec_b[k].elem;
    }

    return ZpMersenneIntElement(v % p);
  }
};

inline std::ostream& operator<<(std::ostream& s, const ZpMersenneIntElement &a) { return s << a.elem; }


class ZpMersenneLongElement {

//private:
public: //TODO return to private after tesing

  static const unsigned long p = 0x1FFFFFFFFFFFFFFF; // 2^61 - 1
  unsigned long elem;

public:

  ZpMersenneLongElement() { elem = 0; }
  ZpMersenneLongElement(unsigned long elem)
  {
    this->elem = elem;
    if (this->elem >= p){

      this->elem = (this->elem & p) + (this->elem >> 61);

      if (this->elem >= p)
        this->elem -= p;
    }
  }
  ZpMersenneLongElement(const ZpMersenneIntElement &rhs) {
    this->elem = rhs.elem;
  }

  inline ZpMersenneLongElement& operator=(const ZpMersenneLongElement &other) {elem = other.elem; return *this; }

  inline bool operator==(const ZpMersenneLongElement &other) const { return (other.elem == elem); }
  inline bool operator!=(const ZpMersenneLongElement &other) const { return !(other.elem == elem); }

  ZpMersenneLongElement operator+(const ZpMersenneLongElement &f2) const
  {
    ZpMersenneLongElement answer;

    answer.elem = (elem + f2.elem);

    if (answer.elem >= p)
      answer.elem -= p;

    return answer;
  }

  ZpMersenneLongElement operator-(const ZpMersenneLongElement &f2) const
  {
    ZpMersenneLongElement answer;

    long temp =  (long)elem - (long)f2.elem;

    if(temp < 0) {
      answer.elem = temp + p;
    }
    else {
      answer.elem = temp;
    }

    return answer;
  }

  ZpMersenneLongElement operator/(const ZpMersenneLongElement &f2) const
  {
    ZpMersenneLongElement answer;

    mpz_t d;
    mpz_t result;
    mpz_t mpz_elem;
    mpz_t mpz_me;
    mpz_init_set_str (d, "2305843009213693951", 10);
    mpz_init(mpz_elem);
    mpz_init(mpz_me);

    mpz_set_ui(mpz_elem, f2.elem);
    mpz_set_ui(mpz_me, elem);

    mpz_init(result);

    mpz_invert ( result, mpz_elem, d );

    mpz_mul (result, result, mpz_me);
    mpz_mod (result, result, d);


    answer.elem = mpz_get_ui(result);

    return answer;
  }

  ZpMersenneLongElement operator*(const ZpMersenneLongElement &f2) const
  {
    ZpMersenneLongElement answer;

    unsigned long long high;
    unsigned long low = _mulx_u64(elem, f2.elem, &high);

    unsigned long low61 = (low & p);
    unsigned long low62to64 = (low >> 61);
    unsigned long highShift3 = (high << 3);

    unsigned long res = low61 + low62to64 + highShift3;

    if (res >= p)
      res -= p;

    answer.elem = res;

    return answer;
  }

  ZpMersenneLongElement& operator+=(const ZpMersenneLongElement &f2)
  {
    elem = (elem + f2.elem);

    if (elem >= p)
      elem -= p;

    return *this;
  }

  ZpMersenneLongElement& operator-=(const ZpMersenneLongElement &f2)
  {
    long temp =  (long)elem - (long)f2.elem;

    if (temp < 0) {
      elem = temp + p;
    }
    else {
      elem = temp;
    }

    return *this;
  }

  ZpMersenneLongElement& operator*=(const ZpMersenneLongElement &f2)
  {
    unsigned long long high;
    unsigned long low = _mulx_u64(elem, f2.elem, &high);

    unsigned long low61 = (low & p);
    unsigned long low61to64 = (low >> 61);
    unsigned long highShift3 = (high << 3);

    unsigned long res = low61 + low61to64 + highShift3;

    if (res >= p)
      res -= p;

    elem = res;

    return *this;
  }

  static ZpMersenneLongElement dotProdct(ZpMersenneLongElement **&mat_a, const std::vector<std::vector<ZpMersenneLongElement>> &mat_b, int row, int col, int length) {
    uint128_t v = 0, t;

    int quotiant = length / 64;
    int remainder = length - (quotiant * 64);

    int qq = 0;
    for (int q = 0; q < quotiant; ++q) { // OK up to very large "length"
      t = 0;
      for (int k = 0; k < 64; ++k) {
        t += (uint128_t)mat_a[row][qq + k].elem * (uint128_t)mat_b[qq + k][col].elem;
      }
      v += t % p;
      qq += 64;
    }

    t = 0;
    for (int k = 0; k < remainder; ++k) {
      t += (uint128_t)mat_a[row][qq + k].elem * (uint128_t)mat_b[qq + k][col].elem;
    }
    v += t % p;

    ZpMersenneLongElement ret;

    ret.elem = v % p;

    return ret;
  }

  static ZpMersenneLongElement dotProdct(ZpMersenneLongElement **&mat_a, const std::vector<ZpMersenneLongElement> &vec_b, int row, int col, int length) {
    uint128_t v = 0, t;

    int quotiant = length / 64;
    int remainder = length - (quotiant * 64);

    int qq = 0;
    for (int q = 0; q < quotiant; ++q) { // OK up to very large "length"
      t = 0;
      for (int k = 0; k < 64; ++k) {
        t += (uint128_t)mat_a[row][qq + k].elem * (uint128_t)vec_b[qq + k].elem;
      }
      v += t % p;
      qq += 64;
    }

    t = 0;
    for (int k = 0; k < remainder; ++k) {
      t += (uint128_t)mat_a[row][qq + k].elem * (uint128_t)vec_b[qq + k].elem;
    }
    v += t % p;

    ZpMersenneLongElement ret;

    ret.elem = v % p;

    return ret;
  }
};

inline std::ostream& operator<<(std::ostream& s, const ZpMersenneLongElement &a) { return s << a.elem; }


#endif // __LZKP_MERSENNE_H_FILE__
