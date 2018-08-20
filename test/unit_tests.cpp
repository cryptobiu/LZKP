
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

  REQUIRE(eq(s.b,st[14].b));
  REQUIRE(eq(s.b, st.getSeed(7).b)); // Test getter
}

TEST_CASE("prover_constructor") {
  Settings s(50, 8, 31, 10, 7, 13);

  Prover p(s);
}

TEST_CASE("verifier_constuctor") {
  auto M = 50;
  Settings s(M, 8, 31, 10, 7, 3);

  Verifier v(s);

  for (auto i = 0; i < M; ++i)
    REQUIRE(v.E_[i] == false);
}

TEST_CASE("prover_r1") {
  Settings s(50, 8, 31, 10, 7, 13);

  Prover p(s);

  p.r1();
}

TEST_CASE("verifier_r2") {
  auto M = 50;
  Settings s(M, 8, 31, 10, 7, 13);

  Verifier v(s);

  block dummy;
  std::vector<bool> e = v.r2(dummy);

  int sum = 0;

  for (auto i = 0; i < M; ++i) {
    if (e[i])
      ++sum;
  }

  REQUIRE(sum == 13);
}

TEST_CASE("prover_r3") {
  auto M = 50, t = 13;
  Settings s(M, 8, 31, 10, 7, t); // M, N, q, m, n, tau

  Prover p(s);
  Verifier v(s);

  block h_gamma = p.r1(); // Run round 1
  std::vector<bool> E = v.r2(h_gamma); // Run round 2

  std::vector<block> seed, omega;
  block h_pi;

  p.r3(E, seed, omega, h_pi); // Run round 3

  REQUIRE(seed.size() == t);
  REQUIRE(omega.size() == M - t);
}

TEST_CASE("verifier_r4") {
  auto M = 50, t = 13, m = 10, n = 10;
  Settings s(M, 8, 31, m, n, t); // M, N, q, m, n, tau

  Prover p(s);
  Verifier v(s);

  block h_gamma = p.r1(); // Run round 1
  std::vector<bool> E = v.r2(h_gamma); // Run round 2

  std::vector<block> seed, omega;
  block h_pi;

  p.r3(E, seed, omega, h_pi); // Run round 3

  std::vector<std::vector<NTL::ZZ_p>> coefficients;
  v.r4(seed, omega, h_pi, coefficients); // Run round 4

  int e_id = 0;
  for (auto e = 0; e < M; ++e) {
    if (E[e]) {
      REQUIRE(eq(p.h_[e].b, v.h_[e].b));
      REQUIRE(coefficients[e_id++].size() == m + n);
    }
  }
}

TEST_CASE("prover_r5") {
  auto M = 50, t = 13, m = 10, n = 10;
  Settings s(M, 8, 31, m, n, t); // M, N, q, m, n, tau

  Prover p(s);
  Verifier v(s);

  block h_gamma = p.r1(); // Run round 1
  std::vector<bool> E = v.r2(h_gamma); // Run round 2

  std::vector<block> seed, omega;
  block h_pi;

  p.r3(E, seed, omega, h_pi); // Run round 3

  std::vector<std::vector<NTL::ZZ_p>> coefficients;
  v.r4(seed, omega, h_pi, coefficients); // Run round 4

  block h_psi = p.r5(coefficients);
}

TEST_CASE("verifier_r6") {
  auto M = 50, t = 13, m = 10, n = 10;
  Settings s(M, 8, 31, m, n, t); // M, N, q, m, n, tau

  Prover p(s);
  Verifier v(s);

  block h_gamma = p.r1(); // Run round 1
  std::vector<bool> E = v.r2(h_gamma); // Run round 2

  std::vector<block> seed, omega;
  block h_pi;

  p.r3(E, seed, omega, h_pi); // Run round 3

  std::vector<std::vector<NTL::ZZ_p>> coefficients;
  v.r4(seed, omega, h_pi, coefficients); // Run round 4

  block h_psi = p.r5(coefficients);

  std::vector<int> i_bar;

  v.r6(h_psi, i_bar);
}
