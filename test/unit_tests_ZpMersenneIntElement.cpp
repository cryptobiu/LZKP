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
#include "sac_prover_logic.h"
#include "sac_verifier_logic.h"
#include "fields/mersenne.h"


using namespace lzkp;


TEST_CASE("cac_prover_logic<ZpMersenneIntElement>_constructor") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, tau, N, n, m);

  CacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);

  REQUIRE(p.M == M);
  REQUIRE(p.tau == tau);
  REQUIRE(p.N == N);
  REQUIRE(p.n == n);
  REQUIRE(p.m == m);
}

TEST_CASE("cac_verifier_logic<ZpMersenneIntElement>_constructor") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, tau, N, n, m);

  CacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

  REQUIRE(v.M == M);
  REQUIRE(v.tau == tau);
  REQUIRE(v.N == N);
  REQUIRE(v.n == n);
  REQUIRE(v.m == m);
}

TEST_CASE("cac_prover_logic<ZpMersenneIntElement>_r1") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, tau, N, n, m);

  CacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);

  block h_gamma;

  p.r1(h_gamma); // Run round 1
}

TEST_CASE("cac_verifier_logic<ZpMersenneIntElement>_r2") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, tau, N, n, m);

  CacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  CacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

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

TEST_CASE("cac_prover_logic<ZpMersenneIntElement>_r3") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, tau, N, n, m);

  CacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  CacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

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

TEST_CASE("cac_verifier_logic<ZpMersenneIntElement>_r4") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, tau, N, n, m);

  CacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  CacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

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

TEST_CASE("cac_prover_logic<ZpMersenneIntElement>_r5") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, tau, N, n, m);

  CacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  CacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

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

TEST_CASE("cac_verifier_logic<ZpMersenneIntElement>_r6") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, tau, N, n, m);

  CacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  CacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

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

TEST_CASE("cac_prover_logic<ZpMersenneIntElement>_r7") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, tau, N, n, m);

  CacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  CacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

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
  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s;
  std::vector<ZpMersenneIntElement> o_i_bar;

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

      if (cur_i_bar >= item.second.first &&  cur_i_bar < item.second.second) { // seed_e_i_bar[e] is a descendant of the current node
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

TEST_CASE("cac_verifier_logic<ZpMersenneIntElement>_r8") {
  auto M = 50, tau = 13, N = 4, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, tau, N, n, m);

  CacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  CacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

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
  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s;
  std::vector<ZpMersenneIntElement> o_i_bar;

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

TEST_CASE("cac_verifier_logic<ZpMersenneIntElement>_r8_false") {
  auto M = 50, tau = 13, N = 4, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  secret[0] += ZpMersenneIntElement(1); // Proof should be rejected now

  Parameters param(M, tau, N, n, m);

  CacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  CacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

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
  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s;
  std::vector<ZpMersenneIntElement> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8

  REQUIRE(!flag);
}

TEST_CASE("sac_prover_logic<ZpMersenneIntElement>_constructor") {
  auto M = 50, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, 0, N, n, m);

  SacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);

  REQUIRE(p.M == M);
  REQUIRE(p.N == N);
  REQUIRE(p.n == n);
  REQUIRE(p.m == m);
}

TEST_CASE("sac_verifier_logic<ZpMersenneIntElement>_constructor") {
  auto M = 50, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, 0, N, n, m);

  SacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

  REQUIRE(v.M == M);
  REQUIRE(v.N == N);
  REQUIRE(v.n == n);
  REQUIRE(v.m == m);
}

TEST_CASE("sac_prover_logic<ZpMersenneIntElement>_r1") {
  auto M = 50, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, 0, N, n, m);

  SacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);

  block h_gamma;

  p.r1(h_gamma); // Run round 1
}

TEST_CASE("sac_verifier_logic<ZpMersenneIntElement>_r2") {
  auto M = 50, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, 0, N, n, m);

  SacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  SacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  block seed_ell;

  v.r2(h_gamma, seed_ell); // Run round 2

  REQUIRE(eq(p.h_gamma_.b, v.h_gamma_.b));
}

TEST_CASE("sac_prover_logic<ZpMersenneIntElement>_r3") {
  auto M = 50, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, 0, N, n, m);

  SacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  SacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  block seed_ell;

  v.r2(h_gamma, seed_ell); // Run round 2

  block h_pi, h_psi, h_theta;

  p.r3(seed_ell, h_pi, h_psi, h_theta); // Run round 3

  REQUIRE(eq(p.seed_ell_.b, v.seed_ell_.b));

  for (auto e = 0; e < M; ++e) {
    REQUIRE(p.provers_[e]->ep_.size() == m);
    REQUIRE(p.provers_[e]->be_.size() == n);
    REQUIRE(p.provers_[e]->ga_.size() == m);
    REQUIRE(p.provers_[e]->de_.size() == m);

    for (auto i = 0; i < m; ++i) {
      REQUIRE(p.provers_[e]->ep_[i] == v.verifiers_[e]->ep_[i]);
      REQUIRE(p.provers_[e]->ga_[i] == v.verifiers_[e]->ga_[i]);
    }

    for (auto i = 0; i < n; ++i) {
      REQUIRE(p.provers_[e]->be_[i] == v.verifiers_[e]->be_[i]);
    }

    ZpMersenneIntElement o_sigma = ZpMersenneIntElement(0);
    for (auto i = 0; i < N; ++i) {
      o_sigma += p.provers_[e]->o_[i];
    }
    REQUIRE(o_sigma == ZpMersenneIntElement(0));
  }
}

TEST_CASE("sac_prover_logic<ZpMersenneIntElement>_r4") {
  auto M = 50, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, 0, N, n, m);

  SacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  SacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  block seed_ell;

  v.r2(h_gamma, seed_ell); // Run round 2

  block h_pi, h_psi, h_theta;

  p.r3(seed_ell, h_pi, h_psi, h_theta); // Run round 3

  std::vector<int> i_bar;

  v.r4(h_pi, h_psi, h_theta, i_bar); // Run round 6

  REQUIRE(eq(p.h_pi_.b, v.h_pi_.b));
  REQUIRE(eq(p.h_psi_.b, v.h_psi_.b));
  REQUIRE(eq(p.h_theta_.b, v.h_theta_.b));

  REQUIRE(i_bar.size() == M);
  for (auto i = 0; i < M; ++i) {
    REQUIRE(i_bar[i] >= 0);
    REQUIRE(i_bar[i] < N);
  }
}

TEST_CASE("sac_prover_logic<ZpMersenneIntElement>_r5") {
  auto M = 50, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, 0, N, n, m);

  SacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  SacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  block seed_ell;

  v.r2(h_gamma, seed_ell); // Run round 2

  block h_pi, h_psi, h_theta;

  p.r3(seed_ell, h_pi, h_psi, h_theta); // Run round 3

  std::vector<int> i_bar;

  v.r4(h_pi, h_psi, h_theta, i_bar); // Run round 4

  block seed_global;
  std::vector<std::vector<block>> seed_tree;
  std::vector<block> gamma_i_bar;
  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s, s_square;
  std::vector<ZpMersenneIntElement> o_i_bar, v_i_bar;

  p.r5(i_bar, seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 5

  for (auto e = 0; e < M; ++e) {
    REQUIRE(p.provers_[e]->i_bar_ == v.verifiers_[e]->i_bar_);
  }

  REQUIRE(eq(seed_global.b, p.seed_global_.b));

  REQUIRE(seed_tree.size() == M);
  for (auto e = 0; e < M; ++e) {
    const auto cur_i_bar = i_bar[e];
    REQUIRE(cur_i_bar == p.provers_[e]->i_bar_);
    REQUIRE(seed_tree[e].size() == (log(N)/log(2)));

    std::stack<std::pair<int, std::pair<int, int>>> seed_stack;

    seed_stack.push(std::make_pair(0, std::make_pair(0, N)));

    int id = 0;
    while (!seed_stack.empty()) {
      auto item = seed_stack.top();
      seed_stack.pop();

      // Check for leaf
      if (item.second.second - item.second.first == 1) {
        if (item.second.first != cur_i_bar)
          REQUIRE(eq(seed_tree[e].at(id++).b, p.provers_[e]->seed_tree_.getSeed(item.second.first).b));
        continue;
      }

      if (cur_i_bar >= item.second.first &&  cur_i_bar < item.second.second) { // seed_e_i_bar[e] is a descendant of the current node
        seed_stack.push(std::make_pair(item.first * 2 + 2, std::make_pair((item.second.first + item.second.second) / 2, item.second.second))); // Add right child
        seed_stack.push(std::make_pair(item.first * 2 + 1, std::make_pair(item.second.first, (item.second.first + item.second.second) / 2))); // Add left child
      }
      else {
        REQUIRE(eq(seed_tree[e].at(id++).b, p.provers_[e]->seed_tree_[item.first].b));
      }
    }

    REQUIRE(eq(gamma_i_bar[e].b, p.provers_[e]->gamma_[cur_i_bar].b));

    REQUIRE(alpha_i_bar[e].size() == m);
    for (auto k = 0; k < m; ++k) {
      REQUIRE(alpha_i_bar[e][k] == p.provers_[e]->alpha_[k][cur_i_bar]);
    }

    REQUIRE(o_i_bar[e] == p.provers_[e]->o_[cur_i_bar]);
    REQUIRE(v_i_bar[e] == p.provers_[e]->v_[cur_i_bar]);

    if (cur_i_bar != N - 1) {
      REQUIRE(b_square[e].size() == m);
      REQUIRE(s[e].size() == m);
      REQUIRE(s_square[e].size() == m);

      for (auto k = 0; k < m; ++k) {
        REQUIRE(b_square[e][k] == p.provers_[e]->b_square_[k][N - 1]);
        REQUIRE(s[e][k] == p.provers_[e]->s_[k][N - 1]);
        REQUIRE(s_square[e][k] == p.provers_[e]->s_square_[k][N - 1]);
      }
    }
  }
}

TEST_CASE("sac_verifier_logic<ZpMersenneIntElement>_r6") {
  auto M = 50, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Parameters param(M, 0, N, n, m);

  SacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  SacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  block seed_ell;

  v.r2(h_gamma, seed_ell); // Run round 2

  block h_pi, h_psi, h_theta;

  p.r3(seed_ell, h_pi, h_psi, h_theta); // Run round 3

  std::vector<int> i_bar;

  v.r4(h_pi, h_psi, h_theta, i_bar); // Run round 4

  block seed_global;
  std::vector<std::vector<block>> seed_tree;
  std::vector<block> gamma_i_bar;
  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s, s_square;
  std::vector<ZpMersenneIntElement> o_i_bar, v_i_bar;

  p.r5(i_bar, seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 5

  bool flag = v.r6(seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 6

  for (auto e = 0; e < M; ++e) {
    REQUIRE(!v.verifiers_[e]->reject_);
  }

  // Check seed reconstruction
  for (auto e = 0; e < M; ++e) {
    const auto cur_i_bar = i_bar[e];

    for (auto i = 0; i < N; ++i) {
      if (i == cur_i_bar)
        continue;

      REQUIRE(eq(p.provers_[e]->seed_tree_.getSeed(i).b, v.verifiers_[e]->partial_seeds_[i].b));

      REQUIRE(eq(p.provers_[e]->r_[i].b, v.verifiers_[e]->r_[i].b));

      for (auto k = 0; k < m; ++k) {
        REQUIRE(p.provers_[e]->b_[k][i] == v.verifiers_[e]->b_[k][i]);

        if (i != N - 1) {
          REQUIRE(p.provers_[e]->b_square_[k][i] == v.verifiers_[e]->b_square_[k][i]);
          REQUIRE(p.provers_[e]->s_[k][i] == v.verifiers_[e]->s_[k][i]);
          REQUIRE(p.provers_[e]->s_square_[k][i] == v.verifiers_[e]->s_square_[k][i]);
        }
      }

      REQUIRE(eq(p.provers_[e]->gamma_[i].b, v.verifiers_[e]->gamma_[i].b));

      REQUIRE(p.provers_[e]->o_[i] == v.verifiers_[e]->o_[i]);

      REQUIRE(eq(p.provers_[e]->seed_global_.b, v.verifiers_[e]->seed_global_.b));

      REQUIRE(p.provers_[e]->v_[i] == v.verifiers_[e]->v_[i]);
    }

    REQUIRE(p.provers_[e]->de_.size() == v.verifiers_[e]->de_.size());
    for (auto k = 0; k < m; ++k) {
      REQUIRE(p.provers_[e]->de_[k] == v.verifiers_[e]->de_[k]);
    }

    REQUIRE(eq(p.provers_[e]->g_.b, v.verifiers_[e]->g_.b));
    REQUIRE(eq(p.provers_[e]->w_.b, v.verifiers_[e]->w_.b));
    REQUIRE(eq(p.provers_[e]->u_.b, v.verifiers_[e]->u_.b));

    REQUIRE(eq(p.provers_[e]->h_.b, v.verifiers_[e]->h_.b));

    REQUIRE(eq(p.provers_[e]->pi_.b, v.verifiers_[e]->pi_.b));

    REQUIRE(eq(p.provers_[e]->psi_.b, v.verifiers_[e]->psi_.b));

    REQUIRE(eq(p.provers_[e]->theta_.b, v.verifiers_[e]->theta_.b));
  }

  REQUIRE(flag);
}

TEST_CASE("sac_verifier_logic<ZpMersenneIntElement>_r6_reject") {
  auto M = 50, N = 8, n = 4, m = 32;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneIntElement>> a;
  std::vector<ZpMersenneIntElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneIntElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneIntElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneIntElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  secret[0] += ZpMersenneIntElement(1); // Proof should be rejected now

  Parameters param(M, 0, N, n, m);

  SacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
  SacVerifierLogic<ZpMersenneIntElement> v(param, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  block seed_ell;

  v.r2(h_gamma, seed_ell); // Run round 2

  block h_pi, h_psi, h_theta;

  p.r3(seed_ell, h_pi, h_psi, h_theta); // Run round 3

  std::vector<int> i_bar;

  v.r4(h_pi, h_psi, h_theta, i_bar); // Run round 4

  block seed_global;
  std::vector<std::vector<block>> seed_tree;
  std::vector<block> gamma_i_bar;
  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s, s_square;
  std::vector<ZpMersenneIntElement> o_i_bar, v_i_bar;

  p.r5(i_bar, seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 5

  bool flag = v.r6(seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 6

  REQUIRE(!flag);
}

//TEST_CASE("sac_verifier_logic<ZpMersenneIntElement>_protocol_big_numbers") {
//  auto M = 700, N = 32, n = 512, m = 4096;
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
//  Parameters param(M, 0, N, n, m);
//
//  SacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
//  SacVerifierLogic<ZpMersenneIntElement> v(param, a, t);
//
//  block h_gamma;
//
//  p.r1(h_gamma); // Run round 1
//
//  block seed_ell;
//
//  v.r2(h_gamma, seed_ell); // Run round 2
//
//  block h_pi, h_psi, h_theta;
//
//  p.r3(seed_ell, h_pi, h_psi, h_theta); // Run round 3
//
//  std::vector<int> i_bar;
//
//  v.r4(h_pi, h_psi, h_theta, i_bar); // Run round 4
//
//  block seed_global;
//  std::vector<std::vector<block>> seed_tree;
//  std::vector<block> gamma_i_bar;
//  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s, s_square;
//  std::vector<ZpMersenneIntElement> o_i_bar, v_i_bar;
//
//  p.r5(i_bar, seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 5
//
//  bool flag = v.r6(seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 6
//
//  REQUIRE(flag);
//}
//
//TEST_CASE("sac_verifier_logic<ZpMersenneIntElement>_protocol_big_numbers_reject") {
//  auto M = 700, N = 32, n = 512, m = 4096;
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
//  secret[0] += ZpMersenneIntElement(1); // Proof should be rejected now
//
//  Parameters param(M, 0, N, n, m);
//
//  SacProverLogic<ZpMersenneIntElement> p(param, a, t, secret);
//  SacVerifierLogic<ZpMersenneIntElement> v(param, a, t);
//
//  block h_gamma;
//
//  p.r1(h_gamma); // Run round 1
//
//  block seed_ell;
//
//  v.r2(h_gamma, seed_ell); // Run round 2
//
//  block h_pi, h_psi, h_theta;
//
//  p.r3(seed_ell, h_pi, h_psi, h_theta); // Run round 3
//
//  std::vector<int> i_bar;
//
//  v.r4(h_pi, h_psi, h_theta, i_bar); // Run round 4
//
//  block seed_global;
//  std::vector<std::vector<block>> seed_tree;
//  std::vector<block> gamma_i_bar;
//  std::vector<std::vector<ZpMersenneIntElement>> alpha_i_bar, b_square, s, s_square;
//  std::vector<ZpMersenneIntElement> o_i_bar, v_i_bar;
//
//  p.r5(i_bar, seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 5
//
//  bool flag = v.r6(seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 6
//
//  REQUIRE(!flag);
//}