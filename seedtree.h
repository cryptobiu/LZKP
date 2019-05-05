#ifndef LZKP_SEEDTREE_H
#define LZKP_SEEDTREE_H

#include <vector>

#include "block.h"
#include <cryptoTools/Crypto/PRNG.h>
#include <algorithm>
#include <functional>


namespace lzkp {


class SeedTree {
public:
    SeedTree(osuCrypto::u32 size = 0);
    SeedTree(osuCrypto::u32 size, const block &master_seed);
    void resize(osuCrypto::u32 size);
    void generate(const block &master_seed);
    block getBlock(int idx);
    block getSeed(int idx); // not using operator []

    const block& operator[](int idx) const;

private:
public: // Makes tests easier
    osuCrypto::u32 size_;
    std::vector<block> seeds_;
    std::vector<osuCrypto::PRNG> prngs_;

    osuCrypto::PRNG prng_;

    bool state;
};


}

#endif //LZKP_SEEDTREE_H
