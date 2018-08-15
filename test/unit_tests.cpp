
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include <seedtree.h>
#include "catch.hpp"

#include <cryptoTools/Crypto/PRNG.h>
#include "seedtree.h"
#include "prover.h"
#include "verifier.h"

using namespace lzkp;

TEST_CASE("demo") {

   REQUIRE(0 == 0);

}

TEST_CASE("seedtree") {
  SeedTree st(8);

  block blk;
  blk.b = osuCrypto::sysRandomSeed();

  st.Generate(blk);

  osuCrypto::PRNG prng;

  block s0 = st[0];

  prng.SetSeed(s0.b);
  s0 = prng.get<block>();

  prng.SetSeed(s0.b);
  s0 = prng.get<block>();

  prng.SetSeed(s0.b);
  s0 = prng.get<block>();

   REQUIRE(eq(s0.b,st[7].b));
  //REQUIRE((s0.halves[0] == st[7].halves[0]));
  //REQUIRE((s0.halves[1] == st[7].halves[1]));
}

TEST_CASE("prover_r1") {
  Settings s(5, 8, 31, 10, 3);

  Prover p(s);

  p.r1();
}

TEST_CASE("verifier_constuctor") {
  Settings s(50, 8, 31, 10, 3);

  Verifier v(s);

//  for (auto i = 0; i < 50; ++i)
//    REQUIRE(v.E_[i] == false);
}

TEST_CASE("verifier_r2") {
  Settings s(50, 8, 31, 10, 13);

  Verifier v(s);

  std::vector<bool> e = v.r2();

  int sum = 0;

  for (auto i = 0; i < 50; ++i) {
    if (e[i])
      ++sum;
  }

  REQUIRE(sum == 13);
}
