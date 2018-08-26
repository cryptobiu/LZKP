
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include <cryptoTools/Crypto/PRNG.h>
#include "seedtree.h"
#include "prover.h"
#include "verifier.h"

#include <stack>
#include <seedtree.h>

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

  block h_psi = p.r5(coefficients); // Run round 5
}

TEST_CASE("verifier_r6") {
  auto M = 50, N = 8, t = 13, m = 64, n = 16;
  Settings s(M, N, 31, m, n, t); // M, N, q, m, n, tau

  Prover p(s);
  Verifier v(s);

  block h_gamma = p.r1(); // Run round 1
  std::vector<bool> E = v.r2(h_gamma); // Run round 2

  std::vector<block> seed, omega;
  block h_pi;

  p.r3(E, seed, omega, h_pi); // Run round 3

  std::vector<std::vector<NTL::ZZ_p>> coefficients;
  v.r4(seed, omega, h_pi, coefficients); // Run round 4

  block h_psi = p.r5(coefficients); // Run round 5

  std::vector<int> i_bar;

  v.r6(h_psi, i_bar); // Run round 6

  REQUIRE(i_bar.size() == M - t);
  for (auto i = 0; i < M - t; ++i) {
    REQUIRE(i_bar[i] >= 0);
    REQUIRE(i_bar[i] < N);
  }
}

TEST_CASE("prover_r7") {
  auto M = 50, N = 8, t = 13, m = 64, n = 16;
  Settings set(M, N, 31, m, n, t); // M, N, q, m, n, tau

  Prover p(set);
  Verifier v(set);

  block h_gamma = p.r1(); // Run round 1
  std::vector<bool> E = v.r2(h_gamma); // Run round 2

  std::vector<block> seed, omega;
  block h_pi;

  p.r3(E, seed, omega, h_pi); // Run round 3

  std::vector<std::vector<NTL::ZZ_p>> coefficients;
  v.r4(seed, omega, h_pi, coefficients); // Run round 4

  block h_psi = p.r5(coefficients); // Run round 5

  std::vector<int> i_bar;

  v.r6(h_psi, i_bar); // Run round 6

  block seed_e_bar;
  std::vector<std::vector<block>> seed_tree;
  std::vector<block> gamma_i_bar;
  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
  std::vector<NTL::ZZ_p> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  REQUIRE(eq(seed_e_bar.b, p.seed_e_bar_.b));

  REQUIRE(seed_tree.size() == M - t);
  int e_it = 0;
  for (auto e = 0; e < M; ++e) {
    if (E[e])
      continue;

    const auto cur_i_bar = i_bar[e_it];
    REQUIRE(seed_tree[e_it].size() == (log(N)/log(2)));

    std::stack<std::pair<int, std::pair<int, int>>> seed_stack;

    seed_stack.push(std::make_pair(0, std::make_pair(0, N)));

    int id = 0;
    while (!seed_stack.empty()) {
      auto item = seed_stack.top();
      seed_stack.pop();

      // Check for leaf
      if (item.second.second - item.second.first == 1) {
        if (item.second.first != cur_i_bar)
          REQUIRE(eq(seed_tree[e_it].at(id++).b, p.seed_tree_[e].getSeed(item.second.first).b));
        continue;
      }

      if (cur_i_bar >= item.second.first &&  cur_i_bar < item.second.second) { // seed_e_i_bar[e] is a descendant of the current nocde
        seed_stack.push(std::make_pair(item.first * 2 + 2, std::make_pair((item.second.first + item.second.second) / 2, item.second.second))); // Add right child
        seed_stack.push(std::make_pair(item.first * 2 + 1, std::make_pair(item.second.first, (item.second.first + item.second.second) / 2))); // Add left child
      }
      else {
        REQUIRE(eq(seed_tree[e_it].at(id++).b, p.seed_tree_[e][item.first].b));
      }
    }

    REQUIRE(eq(gamma_i_bar[e_it].b, p.gamma_[e][cur_i_bar].b));

    REQUIRE(alpha_i_bar[e_it].size() == m);
    for (auto k = 0; k < m; ++k) {
      REQUIRE(alpha_i_bar[e_it][k] == p.alpha_[e][k][cur_i_bar]);
    }

    REQUIRE(o_i_bar[e_it] == p.o_[e_it][cur_i_bar]);

    if (cur_i_bar != N - 1) {
      REQUIRE(b_square[e_it].size() == m);
      REQUIRE(s[e_it].size() == m);

      for (auto k = 0; k < m; ++k) {
        REQUIRE(b_square[e_it][k] == p.b_square_[e][k][N - 1]);
        REQUIRE(s[e_it][k] == p.s_[e][k][N - 1]);
      }
    }

    e_it++;
  }
}

TEST_CASE("verifier_r8") {
  auto M = 50, N = 8, t = 13, m = 64, n = 16;
  Settings set(M, N, 31, m, n, t); // M, N, q, m, n, tau

  Prover p(set);
  Verifier v(set);

  block h_gamma = p.r1(); // Run round 1
  std::vector<bool> E = v.r2(h_gamma); // Run round 2

  std::vector<block> seed, omega;
  block h_pi;

  p.r3(E, seed, omega, h_pi); // Run round 3

  std::vector<std::vector<NTL::ZZ_p>> coefficients;
  v.r4(seed, omega, h_pi, coefficients); // Run round 4

  block h_psi = p.r5(coefficients); // Run round 5

  std::vector<int> i_bar;

  v.r6(h_psi, i_bar); // Run round 6

  block seed_e_bar;
  std::vector<std::vector<block>> seed_tree;
  std::vector<block> gamma_i_bar;
  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
  std::vector<NTL::ZZ_p> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  std::vector<std::vector<block>> partial_seeds;

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s, partial_seeds); // Run round 8

  // Check seed reconstruction
  if (flag)
    REQUIRE(partial_seeds.size() == M);
  /*int e_it = 0;
  for (auto e = 0; e < M; ++e) {
    if (E[e]) {
      REQUIRE(partial_seeds[e].size() == 0);
      continue;
    }

    REQUIRE(partial_seeds[e].size() == N);
    const auto cur_i_bar = i_bar[e_it++];

    for (auto i = 0; i < N; ++i) {
      if (i != cur_i_bar)
        //REQUIRE(eq(p.seed_tree_[e].getSeed(i).b, partial_seeds[e][i].b));
        REQUIRE(p.seed_tree_[e].getSeed(i).halves[0] == partial_seeds[e][i].halves[0]);
    }
  }*/
}
