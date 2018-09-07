//
// Created by lzkp on 9/7/18.
//

#include "catch.hpp"

#include <stack>
#include <cryptoTools/Crypto/PRNG.h>
#include <seedtree.h>
#include "seedtree.h"
#include "cac_prover_logic.h"
#include "cac_verifier_logic.h"
#include "fields/field_15_bit.h"


using namespace lzkp;


TEST_CASE("cac_prover_logic<Field15Bit>_constructor") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Field15Bit>> a;
  std::vector<Field15Bit> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Field15Bit(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Field15Bit(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Field15Bit(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, n, m);

  CacProverLogic<Field15Bit> p(set, a, t, secret);

  REQUIRE(p.M == M);
  REQUIRE(p.tau == tau);
  REQUIRE(p.N == N);
  REQUIRE(p.n == n);
  REQUIRE(p.m == m);
}

TEST_CASE("cac_verifier_logic<Field15Bit>_constructor") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Field15Bit>> a;
  std::vector<Field15Bit> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Field15Bit(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Field15Bit(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Field15Bit(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, n, m);

  CacVerifierLogic<Field15Bit> v(set, a, t);

  REQUIRE(v.M == M);
  REQUIRE(v.tau == tau);
  REQUIRE(v.N == N);
  REQUIRE(v.n == n);
  REQUIRE(v.m == m);
}

TEST_CASE("cac_prover_logic<Field15Bit>_r1") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Field15Bit>> a;
  std::vector<Field15Bit> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Field15Bit(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Field15Bit(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Field15Bit(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, n, m);

  CacProverLogic<Field15Bit> p(set, a, t, secret);

  block h_gamma;

  p.r1(h_gamma); // Run round 1
}

TEST_CASE("cac_verifier_logic<Field15Bit>_r2") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Field15Bit>> a;
  std::vector<Field15Bit> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Field15Bit(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Field15Bit(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Field15Bit(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, n, m);

  CacProverLogic<Field15Bit> p(set, a, t, secret);
  CacVerifierLogic<Field15Bit> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  int sum = 0;

  REQUIRE(E.size() == M);
  for (auto i = 0; i < M; ++i) {
    if (E[i])
      ++sum;
  }

  REQUIRE(sum == tau);
}

TEST_CASE("cac_prover_logic<Field15Bit>_r3") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Field15Bit>> a;
  std::vector<Field15Bit> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Field15Bit(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Field15Bit(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Field15Bit(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, n, m);

  CacProverLogic<Field15Bit> p(set, a, t, secret);
  CacVerifierLogic<Field15Bit> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  std::vector<block> seed, omegaN;
  block h_pi;

  p.r3(E, seed, omegaN, h_pi); // Run round 3

  REQUIRE(seed.size() == tau);
  REQUIRE(omegaN.size() == M - tau);
}

TEST_CASE("cac_verifier_logic<Field15Bit>_r4") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Field15Bit>> a;
  std::vector<Field15Bit> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Field15Bit(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Field15Bit(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Field15Bit(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, n, m);

  CacProverLogic<Field15Bit> p(set, a, t, secret);
  CacVerifierLogic<Field15Bit> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  std::vector<block> seed, omegaN;
  block h_pi;

  p.r3(E, seed, omegaN, h_pi); // Run round 3

  block seed_ell;
  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4

  for (auto e = 0; e < M; ++e) {
    if (E[e]) {
      for (auto i = 0; i < N; ++i) {
        REQUIRE(eq(p.provers_[e]->r_[i].b, v.verifiers_[e]->r_[i].b));

        for (auto mm = 0; mm < m; ++mm) {
          REQUIRE(p.provers_[e]->b_[mm][i] == v.verifiers_[e]->b_[mm][i]);
          REQUIRE(p.provers_[e]->b_square_[mm][i] == v.verifiers_[e]->b_square_[mm][i]);
        }

        REQUIRE(eq(p.provers_[e]->seed_tree_.getSeed(i).b, v.verifiers_[e]->seed_tree_.getSeed(i).b));
        REQUIRE(eq(p.provers_[e]->gamma_[i].b, v.verifiers_[e]->gamma_[i].b));
      }

      REQUIRE(eq(p.provers_[e]->h_.b, v.verifiers_[e]->h_.b));
    }
    else {
      REQUIRE(v.verifiers_[e]->coefficients_.size() == n + m);
    }
  }
}

TEST_CASE("cac_prover_logic<Field15Bit>_r5") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Field15Bit>> a;
  std::vector<Field15Bit> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Field15Bit(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Field15Bit(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Field15Bit(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, n, m);

  CacProverLogic<Field15Bit> p(set, a, t, secret);
  CacVerifierLogic<Field15Bit> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  std::vector<block> seed, omegaN;
  block h_pi;

  p.r3(E, seed, omegaN, h_pi); // Run round 3

  block seed_ell;
  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4

  block h_psi;
  p.r5(seed_ell, h_psi); // Run round 5

  for (auto e = 0; e < M; ++e) {
    if (E[e])
      continue;

    REQUIRE(p.provers_[e]->coefficients_.size() == n + m);

    for (auto i = 0; i < n + m; ++i) {
      REQUIRE(p.provers_[e]->coefficients_[i] == v.verifiers_[e]->coefficients_[i]);
    }
  }
}

TEST_CASE("cac_verifier_logic<Field15Bit>_r6") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Field15Bit>> a;
  std::vector<Field15Bit> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Field15Bit(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Field15Bit(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Field15Bit(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, n, m);

  CacProverLogic<Field15Bit> p(set, a, t, secret);
  CacVerifierLogic<Field15Bit> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  std::vector<block> seed, omegaN;
  block h_pi;

  p.r3(E, seed, omegaN, h_pi); // Run round 3

  block seed_ell;
  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4

  block h_psi;
  p.r5(seed_ell, h_psi); // Run round 5

  std::vector<int> i_bar;

  v.r6(h_psi, i_bar); // Run round 6

  REQUIRE(i_bar.size() == M - tau);
  for (auto i = 0; i < M - tau; ++i) {
    REQUIRE(i_bar[i] >= 0);
    REQUIRE(i_bar[i] < N);
  }

  for (auto e = 0, e_id = 0; e < M; ++e) {
    if (E[e])
      continue;

    REQUIRE(i_bar[e_id++] == v.verifiers_[e]->i_bar_);
  }
}

TEST_CASE("cac_prover_logic<Field15Bit>_r7") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Field15Bit>> a;
  std::vector<Field15Bit> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Field15Bit(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Field15Bit(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Field15Bit(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, n, m);

  CacProverLogic<Field15Bit> p(set, a, t, secret);
  CacVerifierLogic<Field15Bit> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  std::vector<block> seed, omegaN;
  block h_pi;

  p.r3(E, seed, omegaN, h_pi); // Run round 3

  block seed_ell;
  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4

  block h_psi;
  p.r5(seed_ell, h_psi); // Run round 5

  std::vector<int> i_bar;

  v.r6(h_psi, i_bar); // Run round 6

  block seed_e_bar;
  std::vector<std::vector<block>> seed_tree;
  std::vector<block> gamma_i_bar;
  std::vector<std::vector<Field15Bit>> alpha_i_bar, b_square, s;
  std::vector<Field15Bit> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  for (auto e = 0; e < M; ++e) {
    if (E[e])
      continue;

    REQUIRE(p.provers_[e]->i_bar_ == v.verifiers_[e]->i_bar_);
  }

  REQUIRE(eq(seed_e_bar.b, p.seed_e_bar_.b));

  REQUIRE(seed_tree.size() == M - tau);
  for (auto e = 0, e_it = 0; e < M; ++e) {
    if (E[e])
      continue;

    const auto cur_i_bar = i_bar[e_it];
    REQUIRE(cur_i_bar == p.provers_[e]->i_bar_);
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
          REQUIRE(eq(seed_tree[e_it].at(id++).b, p.provers_[e]->seed_tree_.getSeed(item.second.first).b));
        continue;
      }

      if (cur_i_bar >= item.second.first &&  cur_i_bar < item.second.second) { // seed_e_i_bar[e] is a descendant of the current nocde
        seed_stack.push(std::make_pair(item.first * 2 + 2, std::make_pair((item.second.first + item.second.second) / 2, item.second.second))); // Add right child
        seed_stack.push(std::make_pair(item.first * 2 + 1, std::make_pair(item.second.first, (item.second.first + item.second.second) / 2))); // Add left child
      }
      else {
        REQUIRE(eq(seed_tree[e_it].at(id++).b, p.provers_[e]->seed_tree_[item.first].b));
      }
    }

    REQUIRE(eq(gamma_i_bar[e_it].b, p.provers_[e]->gamma_[cur_i_bar].b));

    REQUIRE(alpha_i_bar[e_it].size() == m);
    for (auto k = 0; k < m; ++k) {
      REQUIRE(alpha_i_bar[e_it][k] == p.provers_[e]->alpha_[k][cur_i_bar]);
    }

    REQUIRE(o_i_bar[e_it] == p.provers_[e]->o_[cur_i_bar]);

    if (cur_i_bar != N - 1) {
      REQUIRE(b_square[e_it].size() == m);
      REQUIRE(s[e_it].size() == m);

      for (auto k = 0; k < m; ++k) {
        REQUIRE(b_square[e_it][k] == p.provers_[e]->b_square_[k][N - 1]);
        REQUIRE(s[e_it][k] == p.provers_[e]->s_[k][N - 1]);
      }
    }

    e_it++;
  }
}

TEST_CASE("cac_verifier_logic<Field15Bit>_r8") {
  auto M = 50, tau = 13, N = 4, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Field15Bit>> a;
  std::vector<Field15Bit> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Field15Bit(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Field15Bit(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Field15Bit(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, n, m);

  CacProverLogic<Field15Bit> p(set, a, t, secret);
  CacVerifierLogic<Field15Bit> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  std::vector<block> seed, omegaN;
  block h_pi;

  p.r3(E, seed, omegaN, h_pi); // Run round 3

  block seed_ell;
  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4

  block h_psi;
  p.r5(seed_ell, h_psi); // Run round 5

  std::vector<int> i_bar;

  v.r6(h_psi, i_bar); // Run round 6

  block seed_e_bar;
  std::vector<std::vector<block>> seed_tree;
  std::vector<block> gamma_i_bar;
  std::vector<std::vector<Field15Bit>> alpha_i_bar, b_square, s;
  std::vector<Field15Bit> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8

  for (auto e = 0; e < M; ++e) {
    if (E[e])
      continue;

    REQUIRE(!v.verifiers_[e]->reject_);
  }

  // Check seed reconstruction
  for (auto e = 0, e_it = 0; e < M; ++e) {
    if (E[e])
      continue;

    const auto cur_i_bar = i_bar[e_it];

    for (auto i = 0; i < N; ++i) {
      if (i == cur_i_bar)
        continue;

      REQUIRE(eq(p.provers_[e]->seed_tree_.getSeed(i).b, v.verifiers_[e]->partial_seeds_[i].b));

      REQUIRE(eq(p.provers_[e]->r_[i].b, v.verifiers_[e]->r_[i].b));

      for (auto k = 0; k < m; ++k) {
        REQUIRE(p.provers_[e]->b_[k][i] == v.verifiers_[e]->b_[k][i]);

        if (i != N - 1)
          REQUIRE(p.provers_[e]->b_square_[k][i] == v.verifiers_[e]->b_square_[k][i]);
      }

      REQUIRE(eq(p.provers_[e]->gamma_[i].b, v.verifiers_[e]->gamma_[i].b));
    }

    REQUIRE(eq(p.provers_[e]->h_.b, v.verifiers_[e]->h_.b));

    if (cur_i_bar != N - 1)
      REQUIRE(eq(p.provers_[e]->omegaN_.b, v.verifiers_[e]->omegaN_.b));

    REQUIRE(eq(p.provers_[e]->pi_.b, v.verifiers_[e]->pi_.b));

    REQUIRE(eq(p.provers_[e]->psi_.b, v.verifiers_[e]->psi_.b));

    e_it++;
  }

  REQUIRE(flag);
}

TEST_CASE("cac_verifier_logic<Field15Bit>_r8_reject") {
  auto M = 50, tau = 13, N = 4, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Field15Bit>> a;
  std::vector<Field15Bit> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Field15Bit(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Field15Bit(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Field15Bit(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  secret[0] += Field15Bit(1); // Proof should be rejected now

  Settings set(M, tau, N, n, m);

  CacProverLogic<Field15Bit> p(set, a, t, secret);
  CacVerifierLogic<Field15Bit> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  std::vector<block> seed, omegaN;
  block h_pi;

  p.r3(E, seed, omegaN, h_pi); // Run round 3

  block seed_ell;
  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4

  block h_psi;
  p.r5(seed_ell, h_psi); // Run round 5

  std::vector<int> i_bar;

  v.r6(h_psi, i_bar); // Run round 6

  block seed_e_bar;
  std::vector<std::vector<block>> seed_tree;
  std::vector<block> gamma_i_bar;
  std::vector<std::vector<Field15Bit>> alpha_i_bar, b_square, s;
  std::vector<Field15Bit> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8

  REQUIRE(!flag);
}