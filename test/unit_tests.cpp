
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include <stack>
#include <cryptoTools/Crypto/PRNG.h>
#include <seedtree.h>
#include "seedtree.h"
#include "cac_prover_logic.h"
#include "cac_verifier_logic.h"
#include "fields/mersenne.h"


using namespace lzkp;


TEST_CASE("seedtree_seeds") {
  block blk;
  blk.b = osuCrypto::sysRandomSeed();

  SeedTree st(8, blk); // Tree with 8 leaves

  osuCrypto::PRNG prng;

  block s;

  s  = st[0];
  prng.SetSeed(s.b);
  s = prng.get<block>(); // Internal seed 1

  prng.SetSeed(s.b);
  s = prng.get<block>(); // Internal seed 7

  prng.SetSeed(s.b);
  s = prng.get<block>(); // Internal seed 7

  REQUIRE(eq(s.b,st[7].b));
  REQUIRE(eq(s.b, st.getSeed(0).b)); // Test getter

  s  = st[0];
  prng.SetSeed(s.b);
  prng.get<block>();
  s = prng.get<block>(); // Internal seed 2

  prng.SetSeed(s.b);
  prng.get<block>();
  s = prng.get<block>(); // Internal seed 6

  prng.SetSeed(s.b);
  prng.get<block>();
  s = prng.get<block>(); // Internal seed 14

  REQUIRE(eq(s.b, st[14].b));
  REQUIRE(eq(s.b, st.getSeed(7).b)); // Test getter
}