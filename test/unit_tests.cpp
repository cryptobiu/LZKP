
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







//
//TEST_CASE("verifier_r8_MersenneIntElement") {
//  auto M = 50, tau = 13, N = 4, n = 4, m = 32;
//  uint64_t q = ZpMersenneIntElement::p;
//
////  NTL::ZZ_p::init(NTL::ZZ(q));
//
//  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
//
//  std::vector<std::vector<ZpMersenneIntElement>> a;
//  std::vector<ZpMersenneIntElement> t, secret;
//
//  a.resize(n);
//  for (auto i = 0; i < n; ++i)
//    a[i].resize(m);
//
//  t.resize(n);
//  secret.resize(m);
//
//  // Random matrix A
//  for (auto nn = 0; nn < n; ++nn) {
//    for (auto mm = 0; mm < m; ++mm) {
//      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
//    }
//  }
//
//  // Random vector secret
//  for (auto mm = 0; mm < m; ++mm) {
//    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
//  }
//
//  // Calculate t
//  for (auto nn = 0; nn < n; ++nn) {
//    t[nn] = ZpMersenneIntElement(0);
//
//    for (auto mm = 0; mm < m; ++mm) {
//      t[nn] += a[nn][mm] * secret[mm];
//    }
//  }
//
//  Parameters set(M, tau, N, q, n, m);
//
//  CacProverLogic<ZpMersenneIntElement> p(set, a, t, secret);
//  CacVerifierLogic<ZpMersenneIntElement> v(set, a, t);
//
//  block h_gamma;
//
//  p.r1(h_gamma); // Run round 1
//
//  std::vector<uint8_t> E;
//  v.r2(h_gamma, E); // Run round 2
//
//  for (auto e = 0; e < M; ++e) {
//    if (E[e])
//      continue;
//
////    std::cout << "P E " << e << "\t";
//    for (auto i = 0; i < N; ++i) {
////      std::cout << i << " " << p.provers_[e]->seed_tree_.getSeed(i).halves[0] << std::endl;
//    }
//  }
//  std::vector<block> seed, omegaN;
//  block h_pi;
//
//  p.r3(E, seed, omegaN, h_pi); // Run round 3
//
//  block seed_ell;
//  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
//
//  block h_psi;
//  p.r5(seed_ell, h_psi); // Run round 5
//
//  std::vector<int> i_bar;
//
//  v.r6(h_psi, i_bar); // Run round 6
//
//  block seed_e_bar;
//  std::vector<std::vector<block>> seed_tree;
//  std::vector<block> gamma_i_bar;
//  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s;
//  std::vector<ZpMersenneIntElement> o_i_bar;
//
//  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
//
//  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
//
//  for (auto e = 0; e < M; ++e) {
//    if (E[e])
//      continue;
//
//    REQUIRE(!v.verifiers_[e]->reject_);
//  }
//
//  // Check seed reconstruction
//  for (auto e = 0, e_it = 0; e < M; ++e) {
//    if (E[e])
//      continue;
//
//    const auto cur_i_bar = i_bar[e_it];
//
//    for (auto i = 0; i < N; ++i) {
//      if (i == cur_i_bar)
//        continue;
//
//      REQUIRE(eq(p.provers_[e]->seed_tree_.getSeed(i).b, v.verifiers_[e]->partial_seeds_[i].b));
//
//      REQUIRE(eq(p.provers_[e]->r_[i].b, v.verifiers_[e]->r_[i].b));
//
//      for (auto k = 0; k < m; ++k) {
//        REQUIRE(p.provers_[e]->b_[k][i] == v.verifiers_[e]->b_[k][i]);
//
//        if (i != N - 1)
//          REQUIRE(p.provers_[e]->b_square_[k][i] == v.verifiers_[e]->b_square_[k][i]);
//      }
//
//      REQUIRE(eq(p.provers_[e]->gamma_[i].b, v.verifiers_[e]->gamma_[i].b));
//    }
//
//    REQUIRE(eq(p.provers_[e]->h_.b, v.verifiers_[e]->h_.b));
//
//    if (cur_i_bar != N - 1)
//      REQUIRE(eq(p.provers_[e]->omegaN_.b, v.verifiers_[e]->omegaN_.b));
//
//    REQUIRE(eq(p.provers_[e]->pi_.b, v.verifiers_[e]->pi_.b));
//
//    REQUIRE(eq(p.provers_[e]->psi_.b, v.verifiers_[e]->psi_.b));
//
//    e_it++;
//  }
//  REQUIRE(flag);
//}
//
//TEST_CASE("verifier_r8_MersenneLongElement") {
//  auto M = 50, tau = 13, N = 4, n = 4, m = 32;
//  uint64_t q = ZpMersenneLongElement::p;
//
////  NTL::ZZ_p::init(NTL::ZZ(q));
//
//  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
//
//  std::vector<std::vector<ZpMersenneLongElement>> a;
//  std::vector<ZpMersenneLongElement> t, secret;
//
//  a.resize(n);
//  for (auto i = 0; i < n; ++i)
//    a[i].resize(m);
//
//  t.resize(n);
//  secret.resize(m);
//
//  // Random matrix A
//  for (auto nn = 0; nn < n; ++nn) {
//    for (auto mm = 0; mm < m; ++mm) {
//      a[nn][mm] = ZpMersenneLongElement(prng.get<block>().halves[0]);
//    }
//  }
//
//  // Random vector secret
//  for (auto mm = 0; mm < m; ++mm) {
//    secret[mm] = ZpMersenneLongElement(prng.get<block>().bytes[0] % 2);
//  }
//
//  // Calculate t
//  for (auto nn = 0; nn < n; ++nn) {
//    t[nn] = ZpMersenneLongElement(0);
//
//    for (auto mm = 0; mm < m; ++mm) {
//      t[nn] += a[nn][mm] * secret[mm];
//    }
//  }
//
//  Parameters set(M, tau, N, q, n, m);
//
//  CacProverLogic<ZpMersenneLongElement> p(set, a, t, secret);
//  CacVerifierLogic<ZpMersenneLongElement> v(set, a, t);
//
//  block h_gamma;
//
//  p.r1(h_gamma); // Run round 1
//
//  std::vector<uint8_t> E;
//  v.r2(h_gamma, E); // Run round 2
//
//  for (auto e = 0; e < M; ++e) {
//    if (E[e])
//      continue;
//
////    std::cout << "P E " << e << "\t";
//    for (auto i = 0; i < N; ++i) {
////      std::cout << i << " " << p.provers_[e]->seed_tree_.getSeed(i).halves[0] << std::endl;
//    }
//  }
//  std::vector<block> seed, omegaN;
//  block h_pi;
//
//  p.r3(E, seed, omegaN, h_pi); // Run round 3
//
//  block seed_ell;
//  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
//
//  block h_psi;
//  p.r5(seed_ell, h_psi); // Run round 5
//
//  std::vector<int> i_bar;
//
//  v.r6(h_psi, i_bar); // Run round 6
//
//  block seed_e_bar;
//  std::vector<std::vector<block>> seed_tree;
//  std::vector<block> gamma_i_bar;
//  std::vector<std::vector<ZpMersenneLongElement>> alpha_i_bar, b_square, s;
//  std::vector<ZpMersenneLongElement> o_i_bar;
//
//  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
//
//  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
//
//  for (auto e = 0; e < M; ++e) {
//    if (E[e])
//      continue;
//
//    REQUIRE(!v.verifiers_[e]->reject_);
//  }
//
//  // Check seed reconstruction
//  for (auto e = 0, e_it = 0; e < M; ++e) {
//    if (E[e])
//      continue;
//
//    const auto cur_i_bar = i_bar[e_it];
//
//    for (auto i = 0; i < N; ++i) {
//      if (i == cur_i_bar)
//        continue;
//
//      REQUIRE(eq(p.provers_[e]->seed_tree_.getSeed(i).b, v.verifiers_[e]->partial_seeds_[i].b));
//
//      REQUIRE(eq(p.provers_[e]->r_[i].b, v.verifiers_[e]->r_[i].b));
//
//      for (auto k = 0; k < m; ++k) {
//        REQUIRE(p.provers_[e]->b_[k][i] == v.verifiers_[e]->b_[k][i]);
//
//        if (i != N - 1)
//          REQUIRE(p.provers_[e]->b_square_[k][i] == v.verifiers_[e]->b_square_[k][i]);
//      }
//
//      REQUIRE(eq(p.provers_[e]->gamma_[i].b, v.verifiers_[e]->gamma_[i].b));
//    }
//
//    REQUIRE(eq(p.provers_[e]->h_.b, v.verifiers_[e]->h_.b));
//
//    if (cur_i_bar != N - 1)
//      REQUIRE(eq(p.provers_[e]->omegaN_.b, v.verifiers_[e]->omegaN_.b));
//
//    REQUIRE(eq(p.provers_[e]->pi_.b, v.verifiers_[e]->pi_.b));
//
//    REQUIRE(eq(p.provers_[e]->psi_.b, v.verifiers_[e]->psi_.b));
//
//    e_it++;
//  }
//  REQUIRE(flag);
//}
//
////TEST_CASE("verifier_r8_ringNTL") {
////  auto M = 50, tau = 13, N = 4, n = 4, m = 32;
////  uint64_t q = 2147483647;
////
////  NTL::ZZ_p::init(NTL::ZZ(q));
////
////  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
////
////  std::vector<std::vector<RingNTL>> a;
////  std::vector<RingNTL> t, secret;
////
////  a.resize(n);
////  for (auto i = 0; i < n; ++i)
////    a[i].resize(m);
////
////  t.resize(n);
////  secret.resize(m);
////
////  // Random matrix A
////  for (auto nn = 0; nn < n; ++nn) {
////    for (auto mm = 0; mm < m; ++mm) {
////      a[nn][mm] = RingNTL(prng.get<block>().halves[0]);
////    }
////  }
////
////  // Random vector secret
////  for (auto mm = 0; mm < m; ++mm) {
////    secret[mm] = RingNTL(prng.get<block>().bytes[0] % 2);
////  }
////
////  // Calculate t
////  for (auto nn = 0; nn < n; ++nn) {
////    t[nn] = RingNTL(0);
////
////    for (auto mm = 0; mm < m; ++mm) {
////      t[nn] += a[nn][mm] * secret[mm];
////    }
////  }
////
////  Settings set(M, tau, N, q, n, m);
////
////  ProverLogic<RingNTL> p(set, a, t, secret);
////  VerifierLogic<RingNTL> v(set, a, t);
////
////  block h_gamma;
////
////  p.r1(h_gamma); // Run round 1
////
////  std::vector<uint8_t> E;
////  v.r2(h_gamma, E); // Run round 2
////
////  for (auto e = 0; e < M; ++e) {
////    if (E[e])
////      continue;
////
//////    std::cout << "P E " << e << "\t";
////    for (auto i = 0; i < N; ++i) {
//////      std::cout << i << " " << p.provers_[e]->seed_tree_.getSeed(i).halves[0] << std::endl;
////    }
////  }
////  std::vector<block> seed, omegaN;
////  block h_pi;
////
////  p.r3(E, seed, omegaN, h_pi); // Run round 3
////
////  block seed_ell;
////  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
////
////  block h_psi;
////  p.r5(seed_ell, h_psi); // Run round 5
////
////  std::vector<int> i_bar;
////
////  v.r6(h_psi, i_bar); // Run round 6
////
////  block seed_e_bar;
////  std::vector<std::vector<block>> seed_tree;
////  std::vector<block> gamma_i_bar;
////  std::vector<std::vector<RingNTL>> alpha_i_bar, b_square, s;
////  std::vector<RingNTL> o_i_bar;
////
////  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
////
////  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
////
////  for (auto e = 0; e < M; ++e) {
////    if (E[e])
////      continue;
////
////    REQUIRE(!v.verifiers_[e]->reject_);
////  }
////
////  // Check seed reconstruction
////  for (auto e = 0, e_it = 0; e < M; ++e) {
////    if (E[e])
////      continue;
////
////    const auto cur_i_bar = i_bar[e_it];
////
////    for (auto i = 0; i < N; ++i) {
////      if (i == cur_i_bar)
////        continue;
////
////      REQUIRE(eq(p.provers_[e]->seed_tree_.getSeed(i).b, v.verifiers_[e]->partial_seeds_[i].b));
////
////      REQUIRE(eq(p.provers_[e]->r_[i].b, v.verifiers_[e]->r_[i].b));
////
////      for (auto k = 0; k < m; ++k) {
////        REQUIRE(p.provers_[e]->b_[k][i] == v.verifiers_[e]->b_[k][i]);
////
////        if (i != N - 1)
////          REQUIRE(p.provers_[e]->b_square_[k][i] == v.verifiers_[e]->b_square_[k][i]);
////      }
////
////      REQUIRE(eq(p.provers_[e]->gamma_[i].b, v.verifiers_[e]->gamma_[i].b));
////    }
////
////    REQUIRE(eq(p.provers_[e]->h_.b, v.verifiers_[e]->h_.b));
////
////    if (cur_i_bar != N - 1)
////      REQUIRE(eq(p.provers_[e]->omegaN_.b, v.verifiers_[e]->omegaN_.b));
////
////    REQUIRE(eq(p.provers_[e]->pi_.b, v.verifiers_[e]->pi_.b));
////
////    REQUIRE(eq(p.provers_[e]->psi_.b, v.verifiers_[e]->psi_.b));
////
////    e_it++;
////  }
////  REQUIRE(flag);
////}
//
//TEST_CASE("full_protocol_small_numbers") {
//  auto M = 50, tau = 13, N = 4, n = 4, m = 32;
//  uint64_t q = ZpMersenneIntElement::p;
//
//  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
//
//  std::vector<std::vector<ZpMersenneIntElement>> a;
//  std::vector<ZpMersenneIntElement> t, secret;
//
//  a.resize(n);
//  for (auto i = 0; i < n; ++i)
//    a[i].resize(m);
//
//  t.resize(n);
//  secret.resize(m);
//
//  // Random matrix A
//  for (auto nn = 0; nn < n; ++nn) {
//    for (auto mm = 0; mm < m; ++mm) {
//      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
//    }
//  }
//
//  // Random vector secret
//  for (auto mm = 0; mm < m; ++mm) {
//    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
//  }
//
//  // Calculate t
//  for (auto nn = 0; nn < n; ++nn) {
//    t[nn] = ZpMersenneIntElement(0);
//
//    for (auto mm = 0; mm < m; ++mm) {
//      t[nn] += a[nn][mm] * secret[mm];
//    }
//  }
//
//  Parameters set(M, tau, N, q, n, m);
//
//  CacProverLogic<ZpMersenneIntElement> p(set, a, t, secret);
//  CacVerifierLogic<ZpMersenneIntElement> v(set, a, t);
//
//  block h_gamma;
//
//  p.r1(h_gamma); // Run round 1
//
//  std::vector<uint8_t> E;
//  v.r2(h_gamma, E); // Run round 2
//
//  std::vector<block> seed, omegaN;
//  block h_pi;
//
//  p.r3(E, seed, omegaN, h_pi); // Run round 3
//
//  block seed_ell;
//  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
//
//  block h_psi;
//  p.r5(seed_ell, h_psi); // Run round 5
//
//  std::vector<int> i_bar;
//
//  v.r6(h_psi, i_bar); // Run round 6
//
//  block seed_e_bar;
//  std::vector<std::vector<block>> seed_tree;
//  std::vector<block> gamma_i_bar;
//  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s;
//  std::vector<ZpMersenneIntElement> o_i_bar;
//
//  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
//
//  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
//
//  // Check seed reconstruction
//  REQUIRE(flag == true);
//}
//
//TEST_CASE("full_protocol_small_numbers1") {
//  auto M = 50, tau = 13, N = 8, n = 4, m = 32;
//  uint64_t q = ZpMersenneIntElement::p;
//
//  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
//
//  std::vector<std::vector<ZpMersenneIntElement>> a;
//  std::vector<ZpMersenneIntElement> t, secret;
//
//  a.resize(n);
//  for (auto i = 0; i < n; ++i)
//    a[i].resize(m);
//
//  t.resize(n);
//  secret.resize(m);
//
//  // Random matrix A
//  for (auto nn = 0; nn < n; ++nn) {
//    for (auto mm = 0; mm < m; ++mm) {
//      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
//    }
//  }
//
//  // Random vector secret
//  for (auto mm = 0; mm < m; ++mm) {
//    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
//  }
//
//  // Calculate t
//  for (auto nn = 0; nn < n; ++nn) {
//    t[nn] = ZpMersenneIntElement(0);
//
//    for (auto mm = 0; mm < m; ++mm) {
//      t[nn] += a[nn][mm] * secret[mm];
//    }
//  }
//
//  Parameters set(M, tau, N, q, n, m);
//
//  CacProverLogic<ZpMersenneIntElement> p(set, a, t, secret);
//  CacVerifierLogic<ZpMersenneIntElement> v(set, a, t);
//
//  block h_gamma;
//
//  p.r1(h_gamma); // Run round 1
//
//  std::vector<uint8_t> E;
//  v.r2(h_gamma, E); // Run round 2
//
//  std::vector<block> seed, omegaN;
//  block h_pi;
//
//  p.r3(E, seed, omegaN, h_pi); // Run round 3
//
//  block seed_ell;
//  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
//
//  block h_psi;
//  p.r5(seed_ell, h_psi); // Run round 5
//
//  std::vector<int> i_bar;
//
//  v.r6(h_psi, i_bar); // Run round 6
//
//  block seed_e_bar;
//  std::vector<std::vector<block>> seed_tree;
//  std::vector<block> gamma_i_bar;
//  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s;
//  std::vector<ZpMersenneIntElement> o_i_bar;
//
//  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
//
//  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
//
//  // Check seed reconstruction
//  REQUIRE(flag == true);
//}
//
//TEST_CASE("full_protocol_small_numbers2") {
//  auto M = 50, tau = 13, N = 8, n = 8, m = 512;
//  uint64_t q = ZpMersenneIntElement::p;
//
//  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
//
//  std::vector<std::vector<ZpMersenneIntElement>> a;
//  std::vector<ZpMersenneIntElement> t, secret;
//
//  a.resize(n);
//  for (auto i = 0; i < n; ++i)
//    a[i].resize(m);
//
//  t.resize(n);
//  secret.resize(m);
//
//  // Random matrix A
//  for (auto nn = 0; nn < n; ++nn) {
//    for (auto mm = 0; mm < m; ++mm) {
//      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
//    }
//  }
//
//  // Random vector secret
//  for (auto mm = 0; mm < m; ++mm) {
//    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
//  }
//
//  // Calculate t
//  for (auto nn = 0; nn < n; ++nn) {
//    t[nn] = ZpMersenneIntElement(0);
//
//    for (auto mm = 0; mm < m; ++mm) {
//      t[nn] += a[nn][mm] * secret[mm];
//    }
//  }
//
//  Parameters set(M, tau, N, q, n, m);
//
//  CacProverLogic<ZpMersenneIntElement> p(set, a, t, secret);
//  CacVerifierLogic<ZpMersenneIntElement> v(set, a, t);
//
//  block h_gamma;
//
//  p.r1(h_gamma); // Run round 1
//
//  std::vector<uint8_t> E;
//  v.r2(h_gamma, E); // Run round 2
//
//  std::vector<block> seed, omegaN;
//  block h_pi;
//
//  p.r3(E, seed, omegaN, h_pi); // Run round 3
//
//  block seed_ell;
//  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
//
//  block h_psi;
//  p.r5(seed_ell, h_psi); // Run round 5
//
//  std::vector<int> i_bar;
//
//  v.r6(h_psi, i_bar); // Run round 6
//
//  block seed_e_bar;
//  std::vector<std::vector<block>> seed_tree;
//  std::vector<block> gamma_i_bar;
//  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s;
//  std::vector<ZpMersenneIntElement> o_i_bar;
//
//  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
//
//  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
//
//  // Check seed reconstruction
//  REQUIRE(flag == true);
//}
//
//TEST_CASE("full_protocol_small_numbers2_false") {
//  auto M = 50, tau = 13, N = 8, n = 8, m = 512;
//  uint64_t q = ZpMersenneIntElement::p;
//
//  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
//
//  std::vector<std::vector<ZpMersenneIntElement>> a;
//  std::vector<ZpMersenneIntElement> t, secret;
//
//  a.resize(n);
//  for (auto i = 0; i < n; ++i)
//    a[i].resize(m);
//
//  t.resize(n);
//  secret.resize(m);
//
//  // Random matrix A
//  for (auto nn = 0; nn < n; ++nn) {
//    for (auto mm = 0; mm < m; ++mm) {
//      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
//    }
//  }
//
//  // Random vector secret
//  for (auto mm = 0; mm < m; ++mm) {
//    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
//  }
//
//  // Calculate t
//  for (auto nn = 0; nn < n; ++nn) {
//    t[nn] = ZpMersenneIntElement(0);
//
//    for (auto mm = 0; mm < m; ++mm) {
//      t[nn] += a[nn][mm] * secret[mm];
//    }
//  }
//
//  secret[0] += ZpMersenneIntElement(2);
//
//  Parameters set(M, tau, N, q, n, m);
//
//  CacProverLogic<ZpMersenneIntElement> p(set, a, t, secret);
//  CacVerifierLogic<ZpMersenneIntElement> v(set, a, t);
//
//  block h_gamma;
//
//  p.r1(h_gamma); // Run round 1
//
//  std::vector<uint8_t> E;
//  v.r2(h_gamma, E); // Run round 2
//
//  std::vector<block> seed, omegaN;
//  block h_pi;
//
//  p.r3(E, seed, omegaN, h_pi); // Run round 3
//
//  block seed_ell;
//  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
//
//  block h_psi;
//  p.r5(seed_ell, h_psi); // Run round 5
//
//  std::vector<int> i_bar;
//
//  v.r6(h_psi, i_bar); // Run round 6
//
//  block seed_e_bar;
//  std::vector<std::vector<block>> seed_tree;
//  std::vector<block> gamma_i_bar;
//  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s;
//  std::vector<ZpMersenneIntElement> o_i_bar;
//
//  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
//
//  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
//
//  // Check seed reconstruction
//  REQUIRE(flag == false);
//}
//
////TEST_CASE("full_protocol_big_numbers") {
////  auto M = 220, N = 4, tau = 170, m = 4096, n = 1024;
////  uint64_t q = 2305843009213693951L;
////
////  // MOVE THIS LINE TO THE DRIVER...
////  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***
////
////  NTL::Mat<NTL::ZZ_p> a;
////  NTL::Vec<NTL::ZZ_p> t, secret;
////  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
////
////  a.SetDims(n, m); // Fill with values
////  t.SetLength(n); // Fill with values
////  secret.SetLength(m);
////
////  // Random matrix A
////  for (auto nn = 0; nn < n; ++nn) {
////    for (auto mm = 0; mm < m; ++mm) {
////      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
////    }
////  }
////
////  // Random vector secret
////  for (auto mm = 0; mm < m; ++mm) {
////    secret[mm] = prng.get<block>().bytes[0] % 2;
////  }
////
////  // Calculate t
////  for (auto nn = 0; nn < n; ++nn) {
////    t[nn] = 0;
////
////    for (auto mm = 0; mm < m; ++mm) {
////      t[nn] += a[nn][mm] * secret[mm];
////    }
////  }
////
////  Settings set(M, N, q, m, n, tau);
////
////  ProverLogic p(set, a, t, secret);
////  VerifierLogic v(set, a, t);
////
////  block h_gamma;
////
////  p.r1(h_gamma); // Run round 1
////
////  std::vector<bool> E;
////  v.r2(h_gamma, E); // Run round 2
////
////  std::vector<block> seed, omegaN;
////  block h_pi;
////
////  p.r3(E, seed, omegaN, h_pi); // Run round 3
////
////  block seed_ell;
////  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
////
////  block h_psi;
////  p.r5(seed_ell, h_psi); // Run round 5
////
////  std::vector<int> i_bar;
////
////  v.r6(h_psi, i_bar); // Run round 6
////
////  block seed_e_bar;
////  std::vector<std::vector<block>> seed_tree;
////  std::vector<block> gamma_i_bar;
////  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
////  std::vector<NTL::ZZ_p> o_i_bar;
////
////  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
////
////  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
////
////  // Check seed reconstruction
////  REQUIRE(flag == true);
////}
////
////TEST_CASE("full_protocol_big_numbers1") {
////  auto M = 630, N = 64, tau = 615, m = 4096, n = 1024;
////  uint64_t q = 2305843009213693951L;
////
////  // MOVE THIS LINE TO THE DRIVER...
////  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***
////
////  NTL::Mat<NTL::ZZ_p> a;
////  NTL::Vec<NTL::ZZ_p> t, secret;
////  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
////
////  a.SetDims(n, m); // Fill with values
////  t.SetLength(n); // Fill with values
////  secret.SetLength(m);
////
////  // Random matrix A
////  for (auto nn = 0; nn < n; ++nn) {
////    for (auto mm = 0; mm < m; ++mm) {
////      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
////    }
////  }
////
////  // Random vector secret
////  for (auto mm = 0; mm < m; ++mm) {
////    secret[mm] = prng.get<block>().bytes[0] % 2;
////  }
////
////  // Calculate t
////  for (auto nn = 0; nn < n; ++nn) {
////    t[nn] = 0;
////
////    for (auto mm = 0; mm < m; ++mm) {
////      t[nn] += a[nn][mm] * secret[mm];
////    }
////  }
////
////  Settings set(M, N, q, m, n, tau);
////
////  ProverLogic p(set, a, t, secret);
////  VerifierLogic v(set, a, t);
////
////  block h_gamma;
////
////  p.r1(h_gamma); // Run round 1
////
////  std::vector<bool> E;
////  v.r2(h_gamma, E); // Run round 2
////
////  std::vector<block> seed, omegaN;
////  block h_pi;
////
////  p.r3(E, seed, omegaN, h_pi); // Run round 3
////
////  block seed_ell;
////  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
////
////  block h_psi;
////  p.r5(seed_ell, h_psi); // Run round 5
////
////  std::vector<int> i_bar;
////
////  v.r6(h_psi, i_bar); // Run round 6
////
////  block seed_e_bar;
////  std::vector<std::vector<block>> seed_tree;
////  std::vector<block> gamma_i_bar;
////  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
////  std::vector<NTL::ZZ_p> o_i_bar;
////
////  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
////
////  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
////
////  // Check seed reconstruction
////  REQUIRE(flag == true);
////}