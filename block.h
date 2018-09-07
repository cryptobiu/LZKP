//
// Created by lzkp on 9/7/18.
//

#ifndef __LZKP_BLOCK_H_FILE__
#define __LZKP_BLOCK_H_FILE__


#include <cryptoTools/Crypto/PRNG.h>


namespace lzkp {


typedef union {
  osuCrypto::u64 halves[2];
  unsigned char bytes[16];
  osuCrypto::block b;
} block;


}


#endif // __LZKP_BLOCK_H_FILE__
