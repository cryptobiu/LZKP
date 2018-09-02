
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include <cryptoTools/Crypto/PRNG.h>
#include "seedtree.h"
#include "prover.h"
#include "verifier.h"

#include <stack>
#include <seedtree.h>

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

  REQUIRE(eq(s.b,st[14].b));
  REQUIRE(eq(s.b, st.getSeed(7).b)); // Test getter
}

TEST_CASE("prover_constructor") {
  auto M = 50, N = 8, tau = 13, m = 32, n = 4;
  uint64_t q = 31;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
}

TEST_CASE("verifier_constuctor") {
  auto M = 50, N = 8, tau = 13, m = 32, n = 4;
  uint64_t q = 31;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Verifier v(set, a, t);
}

TEST_CASE("prover_r1") {
  auto M = 50, N = 8, tau = 13, m = 32, n = 4;
  uint64_t q = 31;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);

  block h_gamma;

  p.r1(h_gamma); // Run round 1
}

TEST_CASE("verifier_r2") {
  auto M = 50, N = 8, tau = 13, m = 32, n = 4;
  uint64_t q = 31;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
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
  auto M = 50, N = 8, tau = 13, m = 32, n = 4;
  uint64_t q = 31;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
  v.r2(h_gamma, E); // Run round 2

  std::vector<block> seed, omegaN;
  block h_pi;

  p.r3(E, seed, omegaN, h_pi); // Run round 3

  REQUIRE(seed.size() == tau);
  REQUIRE(omegaN.size() == M - tau);
}

TEST_CASE("verifier_r4") {
  auto M = 50, N = 8, tau = 13, m = 32, n = 4;
  uint64_t q = 31;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
  v.r2(h_gamma, E); // Run round 2

  std::vector<block> seed, omegaN;
  block h_pi;

  p.r3(E, seed, omegaN, h_pi); // Run round 3

  //std::vector<std::vector<NTL::ZZ_p>> coefficients;
  block seed_ell;
  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4

  REQUIRE(v.coefficients_.size() == M - tau);
  int c_id = 0;
  for (auto e = 0; e < M; ++e) {
    if (E[e]) {
      REQUIRE(eq(p.h_[e].b, v.h_[e].b));
    }
    else {
      REQUIRE(v.coefficients_[c_id].size() == n + m);
      c_id++;
    }
  }
}

TEST_CASE("prover_r5") {
  auto M = 50, N = 8, tau = 13, m = 32, n = 4;
  uint64_t q = 31;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
  v.r2(h_gamma, E); // Run round 2

  std::vector<block> seed, omegaN;
  block h_pi;

  p.r3(E, seed, omegaN, h_pi); // Run round 3

  block seed_ell;
  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4

  block h_psi;
  p.r5(seed_ell, h_psi); // Run round 5

  REQUIRE(p.coefficients_.size() == M - tau);
  for (auto e = 0; e < M - tau; ++e) {
    REQUIRE(p.coefficients_[e].size() == n + m);

    for (auto i = 0; i < n + m; ++i) {
      REQUIRE(p.coefficients_[e][i] == v.coefficients_[e][i]);
    }
  }
}

TEST_CASE("verifier_r6") {
  auto M = 50, N = 8, tau = 13, m = 32, n = 4;
  uint64_t q = 31;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
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
}

TEST_CASE("prover_r7") {
  auto M = 50, N = 8, tau = 13, m = 32, n = 4;
  uint64_t q = 31;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
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
  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
  std::vector<NTL::ZZ_p> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  REQUIRE(eq(seed_e_bar.b, p.seed_e_bar_.b));

  REQUIRE(seed_tree.size() == M - tau);
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
  auto M = 50, N = 8, tau = 13, m = 32, n = 4;
  uint64_t q = 31;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
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

TEST_CASE("full_protocol_small_numbers") {
  auto M = 50, N = 8, tau = 13, m = 512, n = 32;
  uint64_t q = 32719;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
//  std::cout << "A" << std::endl;
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
//      std::cout << a[nn][mm] << "   ";
    }
//    std::cout << std::endl;
  }

  // Random vector secret
//  std::cout << "Secret" << std::endl;
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
//    std::cout << secret[mm] << "   ";
  }
//  std::cout << std::endl;

  // Calculate t
//  std::cout << "t" << std::endl;
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
//    std::cout << t[nn] << "   ";
  }
//  std::cout << std::endl;

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
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
  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
  std::vector<NTL::ZZ_p> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  std::vector<std::vector<block>> partial_seeds;

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s, partial_seeds); // Run round 8

  // Check seed reconstruction
  REQUIRE(flag == true);
}

TEST_CASE("full_protocol_small_numbers1") {
  auto M = 50, N = 8, tau = 13, m = 512, n = 32;
  uint64_t q = 1048517;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
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
  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
  std::vector<NTL::ZZ_p> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  std::vector<std::vector<block>> partial_seeds;

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s, partial_seeds); // Run round 8

  // Check seed reconstruction
  REQUIRE(flag == true);
}

TEST_CASE("full_protocol_small_numbers2") {
  auto M = 50, N = 8, tau = 13, m = 512, n = 32;
  uint64_t q = 2305843009213693951L;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
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
  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
  std::vector<NTL::ZZ_p> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  std::vector<std::vector<block>> partial_seeds;

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s, partial_seeds); // Run round 8

  // Check seed reconstruction
  REQUIRE(flag == true);
}

TEST_CASE("full_protocol_big_numbers") {
  auto M = 220, N = 4, tau = 170, m = 4096, n = 1024;
  uint64_t q = 2305843009213693951L;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
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
  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
  std::vector<NTL::ZZ_p> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  std::vector<std::vector<block>> partial_seeds;

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s, partial_seeds); // Run round 8

  // Check seed reconstruction
  REQUIRE(flag == true);
}

TEST_CASE("full_protocol_big_numbers1") {
  auto M = 630, N = 64, tau = 615, m = 4096, n = 1024;
  uint64_t q = 2305843009213693951L;

  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  NTL::Mat<NTL::ZZ_p> a;
  NTL::Vec<NTL::ZZ_p> t, secret;
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a.SetDims(n, m); // Fill with values
  t.SetLength(n); // Fill with values
  secret.SetLength(m);

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
      a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
    }
  }

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret[mm] = prng.get<block>().bytes[0] % 2;
  }

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t[nn] = 0;

    for (auto mm = 0; mm < m; ++mm) {
      t[nn] += a[nn][mm] * secret[mm];
    }
  }

  Settings set(M, N, q, m, n, tau);

  Prover p(set, a, t, secret);
  Verifier v(set, a, t);

  block h_gamma;

  p.r1(h_gamma); // Run round 1

  std::vector<bool> E;
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
  std::vector<std::vector<NTL::ZZ_p>> alpha_i_bar, b_square, s;
  std::vector<NTL::ZZ_p> o_i_bar;

  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7

  std::vector<std::vector<block>> partial_seeds;

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s, partial_seeds); // Run round 8

  // Check seed reconstruction
  REQUIRE(flag == true);
}