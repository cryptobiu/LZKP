//
// Created by lzkp on 9/13/18.
//

#ifndef __LZKP_FIELD_59_BIT_H_FILE__
#define __LZKP_FIELD_59_BIT_H_FILE__


#include <iostream>
#include <cstdint>
#include <x86intrin.h>
#include <gmp.h>


typedef unsigned __int128 uint128_t;

struct Field59Bit {
//private:
public: //TODO return to private after tesing
  static const uint64_t p = 0x7FFFFFFFFFFFFC9; // 2^59 - 55
  uint64_t elem;

public:
  Field59Bit() { elem = 0; }
  Field59Bit(uint64_t elem) {
    this->elem = elem % p;
  }

  Field59Bit& operator=(const Field59Bit &other) { elem = other.elem; return *this; }
  bool operator==(const Field59Bit &other) const { return (other.elem == elem); }
  bool operator!=(const Field59Bit &other) const { return !(other.elem == elem); }

  Field59Bit operator+(const Field59Bit &f2) const {
    uint64_t temp = (uint64_t)elem + (uint64_t)f2.elem;

    if (temp >= p) {
      temp -= p;
    }

    return Field59Bit(temp);
  }

  Field59Bit operator-(const Field59Bit &f2) const {
    int64_t temp = (int64_t)elem - (int64_t)f2.elem;

    if (temp < 0) {
      return Field59Bit(temp + p);
    }
    else {
      return Field59Bit(temp);
    }
  }

  Field59Bit operator/(const Field59Bit &f2) const {
      Field59Bit answer;

      mpz_t d;
      mpz_t result;
      mpz_t mpz_elem;
      mpz_t mpz_me;
      mpz_init_set_str (d, "7FFFFFFFFFFFFC9", 16);
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

  Field59Bit operator*(const Field59Bit &f2) const {
    return Field59Bit(((uint128_t)elem * (uint128_t)f2.elem) % p);
  }

  Field59Bit& operator+=(const Field59Bit &f2) {
    uint64_t temp = (uint64_t)elem + (uint64_t)f2.elem;

    if (temp >= p) {
      temp -= p;
    }

    elem = temp;

    return *this;
  }

  Field59Bit& operator-=(const Field59Bit &f2) {
    int64_t temp = (int64_t)elem - (int64_t)f2.elem;

    if (temp < 0) {
      elem = temp + p;
    }
    else {
      elem = temp;
    }

    return *this;
  }

  Field59Bit& operator*=(const Field59Bit &f2)
  {
    elem = ((uint128_t)elem * (uint128_t)f2.elem) % p;

    return *this;
  }
};


inline std::ostream& operator<<(std::ostream& s, const Field59Bit &a) { return s << a.elem; }


#endif // __LZKP_FIELD_59_BIT_H_FILE__
