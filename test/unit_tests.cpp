
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include <stack>
#include <cryptoTools/Crypto/PRNG.h>
#include <seedtree.h>
#include "seedtree.h"
#include "prover_wrapper.h"
#include "verifier_wrapper.h"
#include "Mersenne.h"
#include "Ring31.h"
#include "RingNTL.h"


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

TEST_CASE("prover_constructor") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
}

TEST_CASE("verifier_constuctor") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);
}

TEST_CASE("prover_r1") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);

  block h_gamma;

  p.r1(h_gamma); // Run round 1
}

TEST_CASE("verifier_r2") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);

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

TEST_CASE("prover_r3") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);

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
//
TEST_CASE("verifier_r4") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);

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

TEST_CASE("prover_r5") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);

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

TEST_CASE("verifier_r6") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);

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

TEST_CASE("prover_r7") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);

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

TEST_CASE("verifier_r8") {
  auto M = 2, tau = 1, N = 4, n = 4, m = 32;
  uint64_t q = ZpMersenneLongElement::p;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneLongElement>> a;
  std::vector<ZpMersenneLongElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneLongElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneLongElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneLongElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneLongElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneLongElement> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  for (auto e = 0; e < M; ++e) {
    if (E[e])
      continue;

//    std::cout << "P E " << e << "\t";
    for (auto i = 0; i < N; ++i) {
//      std::cout << i << " " << p.provers_[e]->seed_tree_.getSeed(i).halves[0] << std::endl;
    }
  }
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
  std::vector<std::vector<ZpMersenneLongElement>> alpha_i_bar, b_square, s;
  std::vector<ZpMersenneLongElement> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8

  for (auto e = 0; e < M; ++e) {
    if (E[e])
      continue;

//    REQUIRE(!v.verifiers_[e]->reject_);
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

TEST_CASE("verifier_r8_ring31") {
  auto M = 2, tau = 1, N = 4, n = 4, m = 32;
  uint64_t q = Ring31::p;

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<Ring31>> a;
  std::vector<Ring31> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = Ring31(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = Ring31(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = Ring31(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<Ring31> p(set, a, t, secret);
  VerifierWrapper<Ring31> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  for (auto e = 0; e < M; ++e) {
    if (E[e])
      continue;

//    std::cout << "P E " << e << "\t";
    for (auto i = 0; i < N; ++i) {
//      std::cout << i << " " << p.provers_[e]->seed_tree_.getSeed(i).halves[0] << std::endl;
    }
  }
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
  std::vector<std::vector<Ring31>> alpha_i_bar, b_square, s;
  std::vector<Ring31> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8

  for (auto e = 0; e < M; ++e) {
    if (E[e])
      continue;

//    REQUIRE(!v.verifiers_[e]->reject_);
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

TEST_CASE("verifier_r8_MersenneIntElement") {
  auto M = 50, tau = 13, N = 4, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

//  NTL::ZZ_p::init(NTL::ZZ(q));

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  for (auto e = 0; e < M; ++e) {
    if (E[e])
      continue;

//    std::cout << "P E " << e << "\t";
    for (auto i = 0; i < N; ++i) {
//      std::cout << i << " " << p.provers_[e]->seed_tree_.getSeed(i).halves[0] << std::endl;
    }
  }
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

TEST_CASE("verifier_r8_MersenneLongElement") {
  auto M = 50, tau = 13, N = 4, n = 4, m = 32;
  uint64_t q = ZpMersenneLongElement::p;

//  NTL::ZZ_p::init(NTL::ZZ(q));

  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  std::vector<std::vector<ZpMersenneLongElement>> a;
  std::vector<ZpMersenneLongElement> t, secret;

  a.resize(n);
  for (auto i = 0; i < n; ++i)
    a[i].resize(m);

  t.resize(n);
  secret.resize(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = ZpMersenneLongElement(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = ZpMersenneLongElement(prng.get<block>().bytes[0] % 2);
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = ZpMersenneLongElement(0);

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneLongElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneLongElement> v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2

  for (auto e = 0; e < M; ++e) {
    if (E[e])
      continue;

//    std::cout << "P E " << e << "\t";
    for (auto i = 0; i < N; ++i) {
//      std::cout << i << " " << p.provers_[e]->seed_tree_.getSeed(i).halves[0] << std::endl;
    }
  }
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
  std::vector<std::vector<ZpMersenneLongElement>> alpha_i_bar, b_square, s;
  std::vector<ZpMersenneLongElement> o_i_bar;

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

//TEST_CASE("verifier_r8_ringNTL") {
//  auto M = 50, tau = 13, N = 4, n = 4, m = 32;
//  uint64_t q = 2147483647;
//
//  NTL::ZZ_p::init(NTL::ZZ(q));
//
//  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
//
//  std::vector<std::vector<RingNTL>> a;
//  std::vector<RingNTL> t, secret;
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
//      a[nn][mm] = RingNTL(prng.get<block>().halves[0]);
//    }
//  }
//
//  // Random vector secret
//  for (auto mm = 0; mm < m; ++mm) {
//    secret[mm] = RingNTL(prng.get<block>().bytes[0] % 2);
//  }
//
//  // Calculate t
//  for (auto nn = 0; nn < n; ++nn) {
//    t[nn] = RingNTL(0);
//
//    for (auto mm = 0; mm < m; ++mm) {
//      t[nn] += a[nn][mm] * secret[mm];
//    }
//  }
//
//  Settings set(M, tau, N, q, n, m);
//
//  ProverWrapper<RingNTL> p(set, a, t, secret);
//  VerifierWrapper<RingNTL> v(set, a, t);
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
//  std::vector<std::vector<RingNTL>> alpha_i_bar, b_square, s;
//  std::vector<RingNTL> o_i_bar;
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

TEST_CASE("full_protocol_small_numbers") {
  auto M = 50, tau = 13, N = 4, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);

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

  // Check seed reconstruction
  REQUIRE(flag == true);
}

TEST_CASE("full_protocol_small_numbers1") {
  auto M = 50, tau = 13, N = 8, n = 4, m = 32;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);

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

  // Check seed reconstruction
  REQUIRE(flag == true);
}

TEST_CASE("full_protocol_small_numbers2") {
  auto M = 50, tau = 13, N = 8, n = 8, m = 512;
  uint64_t q = ZpMersenneIntElement::p;

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

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);

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

  // Check seed reconstruction
  REQUIRE(flag == true);
}

TEST_CASE("full_protocol_small_numbers2_false") {
  auto M = 50, tau = 13, N = 8, n = 8, m = 512;
  uint64_t q = ZpMersenneIntElement::p;

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

  secret[0] += ZpMersenneIntElement(2);

  Settings set(M, tau, N, q, n, m);

  ProverWrapper<ZpMersenneIntElement> p(set, a, t, secret);
  VerifierWrapper<ZpMersenneIntElement> v(set, a, t);

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

  // Check seed reconstruction
  REQUIRE(flag == false);
}

//TEST_CASE("full_protocol_big_numbers") {
//  auto M = 220, N = 4, tau = 170, m = 4096, n = 1024;
//  uint64_t q = 2305843009213693951L;
//
//  // MOVE THIS LINE TO THE DRIVER...
//  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***
//
//  NTL::Mat<NTL::ZZ_p> a;
//  NTL::Vec<NTL::ZZ_p> t, secret;
//  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
//
//  a.SetDims(n, m); // Fill with values
//  t.SetLength(n); // Fill with values
//  secret.SetLength(m);
//
//  // Random matrix A
//  for (auto nn = 0; nn < n; ++nn) {
//    for (auto mm = 0; mm < m; ++mm) {
//      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
//    }
//  }
//
//  // Random vector secret
//  for (auto mm = 0; mm < m; ++mm) {
//    secret[mm] = prng.get<block>().bytes[0] % 2;
//  }
//
//  // Calculate t
//  for (auto nn = 0; nn < n; ++nn) {
//    t[nn] = 0;
//
//    for (auto mm = 0; mm < m; ++mm) {
//      t[nn] += a[nn][mm] * secret[mm];
//    }
//  }
//
//  Settings set(M, N, q, m, n, tau);
//
//  ProverWrapper p(set, a, t, secret);
//  VerifierWrapper v(set, a, t);
//
//  block h_gamma;
//
//  p.r1(h_gamma); // Run round 1
//
//  std::vector<bool> E;
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
//  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
//  std::vector<NTL::ZZ_p> o_i_bar;
//
//  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
//
//  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
//
//  // Check seed reconstruction
//  REQUIRE(flag == true);
//}
//
//TEST_CASE("full_protocol_big_numbers1") {
//  auto M = 630, N = 64, tau = 615, m = 4096, n = 1024;
//  uint64_t q = 2305843009213693951L;
//
//  // MOVE THIS LINE TO THE DRIVER...
//  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***
//
//  NTL::Mat<NTL::ZZ_p> a;
//  NTL::Vec<NTL::ZZ_p> t, secret;
//  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());
//
//  a.SetDims(n, m); // Fill with values
//  t.SetLength(n); // Fill with values
//  secret.SetLength(m);
//
//  // Random matrix A
//  for (auto nn = 0; nn < n; ++nn) {
//    for (auto mm = 0; mm < m; ++mm) {
//      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
//    }
//  }
//
//  // Random vector secret
//  for (auto mm = 0; mm < m; ++mm) {
//    secret[mm] = prng.get<block>().bytes[0] % 2;
//  }
//
//  // Calculate t
//  for (auto nn = 0; nn < n; ++nn) {
//    t[nn] = 0;
//
//    for (auto mm = 0; mm < m; ++mm) {
//      t[nn] += a[nn][mm] * secret[mm];
//    }
//  }
//
//  Settings set(M, N, q, m, n, tau);
//
//  ProverWrapper p(set, a, t, secret);
//  VerifierWrapper v(set, a, t);
//
//  block h_gamma;
//
//  p.r1(h_gamma); // Run round 1
//
//  std::vector<bool> E;
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
//  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
//  std::vector<NTL::ZZ_p> o_i_bar;
//
//  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
//
//  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
//
//  // Check seed reconstruction
//  REQUIRE(flag == true);
//}