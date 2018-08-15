//
// Created by roee on 7/26/18.
//

#include "seedtree.h"

using namespace lzkp;

SeedTree::SeedTree(osuCrypto::u32 size) : size_(size), state(0) {
}

void SeedTree::Resize(osuCrypto::u32 size) {
  size_ = size;

  state = 0; // Invalid
}

void SeedTree::Generate(block master_seed) {
  seeds_.resize(2 * size_ - 1);

  prng_.SetSeed(master_seed.b);

  seeds_[0].b = prng_.getSeed();

  for (auto i = 0; i < size_ - 2; ++i) {
    prng_.SetSeed(seeds_[i].b);

    seeds_[2 * i + 1] = prng_.get<block>();
    seeds_[2 * i + 2] = prng_.get<block>();
  }

  state = 1; // Valid
}

const block& SeedTree::operator[](int idx) const {
  return seeds_[idx];
}
