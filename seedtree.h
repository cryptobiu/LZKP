//
// Created by roee on 7/26/18.
//

#ifndef LZKP_SEEDTREE_H
#define LZKP_SEEDTREE_H

#include <vector>

#include <cryptoTools/Crypto/PRNG.h>

namespace lzkp {

typedef union {
    osuCrypto::u64 halves[2];
    unsigned char bytes[16];
    osuCrypto::block b;
} block;

class SeedTree {
public:
    SeedTree(osuCrypto::u32 size = 0);
    void Resize(osuCrypto::u32 size);
    void Generate(block master_seed);

    const block& operator[](int idx) const;

private:
    osuCrypto::u32 size_;
    std::vector<block> seeds_;

    osuCrypto::PRNG prng_;

    bool state;
};
}

#endif //LZKP_SEEDTREE_H
