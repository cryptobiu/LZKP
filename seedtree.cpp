#include "seedtree.h"


using namespace lzkp;

SeedTree::SeedTree(osuCrypto::u32 size) : size_(size), state(0) {
}

SeedTree::SeedTree(osuCrypto::u32 size, const block &master_seed) : size_(size) {
  generate(master_seed);
}

void SeedTree::resize(osuCrypto::u32 size) {
  size_ = size;

  state = 0; // Invalid
}

void SeedTree::generate(const block &master_seed) {
  seeds_.resize(2 * size_ - 1);

  prng_.SetSeed(master_seed.b);

  seeds_[0].b = prng_.getSeed();

  for (auto i = 0u; i < size_ - 1; ++i) {
    prng_.SetSeed(seeds_[i].b);

    seeds_[2 * i + 1] = prng_.get<block>();
    seeds_[2 * i + 2] = prng_.get<block>();
  }

  prngs_.resize(size_);
  for (auto i = 0u; i < size_; ++i) {
    prngs_[i].SetSeed(seeds_[size_ - 1 + i].b);
  }

  state = 1; // Valid
}

block SeedTree::getBlock(int idx) {
  if (idx < 0 || (unsigned )idx >= size_)
    exit(1);

  return prngs_[idx].get<block>();
}

block SeedTree::getSeed(int idx) {
  if (idx < 0 || (unsigned )idx >= size_)
    exit(1);

  return seeds_[size_ - 1 + idx];
}

const block& SeedTree::operator[](int idx) const {
  return seeds_[idx];
}

