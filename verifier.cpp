//
// Created by roee on 7/26/18.
//

#include "verifier.h"

#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Crypto/sha1.h>

#include "seedtree.h"

#include <stack>
#include <queue>

using namespace lzkp;

Verifier::Verifier(const Settings &s) : M(s.M), N(s.N), q(s.q), m(s.m), n(s.n), tau(s.tau) {
  E_.resize(M);

  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  a_.SetDims(n, m); // Fill with values
  t_.SetLength(m); // Fill with values
}

const std::vector<bool> &Verifier::r2(const block &h_gamma) {
  h_gamma_ = h_gamma;

  osuCrypto::PRNG prng;

  prng.SetSeed(osuCrypto::sysRandomSeed());

  for (auto i = M - tau + 1; i <= M; ++i) {
    int r = prng.get<osuCrypto::u32>() % i;

    if (!E_[r])
      E_[r] = true;
    else
      E_[i - 1] = true;
  }

  return E_;
}

void Verifier::r4(const std::vector<block> &seed, const std::vector<block> &omega, const block &h_pi, std::vector<std::vector<NTL::ZZ_p>> &coefficients) {
  seed_ = seed;
  omega_ = omega;
  h_pi_ = h_pi;

  coefficients_.resize(M - tau);

  seed_tree_.resize(M);
  r_.resize(M);
  b_.resize(M);
  b_square_.resize(M);
  gamma_.resize(M);
  h_.resize(M);

  int seed_id = 0;
  int coef_id = 0;

  for (auto e = 0; e < M; ++e) {
    if (E_[e]) { // 1
      // *1* - 1.a
      seed_tree_[e].resize(N);
      gamma_[e].resize(N);

      seed_tree_[e].generate(seed[seed_id++]); // Get the relevant seed

      r_[e].resize(N);
      b_[e].resize(m);
      b_square_[e].resize(m);

      for (auto mm = 0; mm < m; ++mm){
        b_[e][mm].resize(N);
        b_square_[e][mm].resize(N);
      }

      // *1* - 1.b
      for (auto i = 0; i < N - 1; ++i) {
        r_[e][i] = seed_tree_[e].getBlock(i);

        for (auto mm = 0; mm < m; ++mm) {
          b_[e][mm][i]        = NTL::ZZ_p(seed_tree_[e].getBlock(i).halves[0]); // NEED TO FIND A WAY TO USE ALL 128 BITS
          b_square_[e][mm][i] = NTL::ZZ_p(seed_tree_[e].getBlock(i).halves[0]);
        }
      }

      // *1* - 1.c
      r_[e][N-1] = seed_tree_[e].getBlock(N - 1);

      for (auto mm = 0; mm < m; ++mm) {
        b_[e][mm][N-1]        = NTL::ZZ_p(seed_tree_[e].getBlock(N - 1).halves[0]);
      }

      // *1* - 1.d
      for (auto k = 0; k < m; ++k) {
        b_square_[e][k][N - 1] = NTL::ZZ_p(0);
        for (auto i = 0; i < N; ++i) {
          b_square_[e][k][N - 1] += b_[e][k][i];
        }
        b_square_[e][k][N - 1] *= b_square_[e][k][N - 1];
        for (auto i = 0; i < N - 1; ++i) {
          b_square_[e][k][N - 1] -= b_square_[e][k][i];
        }
      }

      // *1* - 1.e
      osuCrypto::SHA1 sha_gamma(sizeof(block));
      for (auto i = 0; i < N - 1; ++i) {
        block blk = seed_tree_[e].getSeed(i);
        sha_gamma.Reset();
        sha_gamma.Update(blk);
        sha_gamma.Update(r_[e][i]);
        sha_gamma.Final(gamma_[e][i]);
      }
      block blk = seed_tree_[e].getSeed(N - 1);
      sha_gamma.Reset();
      sha_gamma.Update(blk);
      unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT
      for (auto i = 0; i < m; ++i) {
        NTL::BytesFromZZ(buf, b_square_[e][i][N - 1]._ZZ_p__rep, NTL::NumBytes(b_square_[e][i][N - 1]._ZZ_p__rep));
        sha_gamma.Update(buf, NTL::NumBytes(b_square_[e][i][N - 1]._ZZ_p__rep));
      }
      sha_gamma.Update(r_[e][N - 1]);
      sha_gamma.Final(gamma_[e][N - 1]);

      // *1* - 1.f
      osuCrypto::SHA1 sha_h(sizeof(block));
      for (auto i = 0; i < N; ++i) {
        sha_h.Update(gamma_[e][i]);
      }
      sha_h.Final(h_[e]); // mb need to zero it first
    }
    else { // 2
      osuCrypto::PRNG prng;
      prng.SetSeed(osuCrypto::sysRandomSeed());

      coefficients_[coef_id].resize(n + m);

      for (auto i = 0; i < n; ++i) {
        coefficients_[coef_id][i] = NTL::ZZ_p(prng.get<block>().halves[0]);
      }
      for (auto i = 0; i < m; ++i) {
        coefficients_[coef_id][n + i] = NTL::ZZ_p(prng.get<block>().halves[0]);
      }

      coef_id++;
    }
  }

  coefficients = coefficients_;
}

void Verifier::r6(const block &h_psi, std::vector<int> &i_bar) {
  h_psi_ = h_psi;

  i_bar_.resize(M - tau);

  osuCrypto::PRNG prng;

  prng.SetSeed(osuCrypto::sysRandomSeed());

  for (auto e = 0; e < M - tau; ++e) {
    i_bar_[e] = prng.get<osuCrypto::u32>() % N;
  }

  i_bar = i_bar_;
}

bool Verifier::r8(const block &seed_e_bar, const std::vector<std::vector<block>> &seed_tree, const std::vector<block> &gamma_i_bar, const std::vector<std::vector<NTL::ZZ_p>> &alpha_i_bar,
                  const std::vector<NTL::ZZ_p> &o_i_bar, const std::vector<std::vector<NTL::ZZ_p>> &b_square,
const std::vector<std::vector<NTL::ZZ_p>> &s, std::vector<std::vector<block>> &partial_seeds) {
  seed_e_bar_ = seed_e_bar;
  prng_e_bar_.SetSeed(seed_e_bar_.b);

  int e_it = 0;

  partial_seeds_.resize(M);
  g_.resize(M - tau);
  w_.resize(M - tau);
  psi_.resize(M - tau);
  o_.resize(M - tau);

  osuCrypto::SHA1 sha_h_pi(sizeof(block));
  osuCrypto::SHA1 sha_h_psi(sizeof(block));

  for (auto e = 0; e < M; ++e) {
    if (E_[e])
      continue;

    // 1.a
    partial_seeds_[e].resize(N);
    int s_it = 0;

    const auto cur_i_bar = i_bar_[e_it];

    std::stack<std::pair<int, std::pair<int, int>>> seed_stack;

    seed_stack.push(std::make_pair(0, std::make_pair(0, N)));

    while (!seed_stack.empty()) {
      auto item = seed_stack.top();
      seed_stack.pop();

      // Check for leaf
      if (item.second.second - item.second.first == 1) {
        if (item.second.first != cur_i_bar) {
          partial_seeds_[e][item.second.first] = seed_tree[e_it][s_it++];
        }
        continue;
      }

      if (cur_i_bar >= item.second.first && cur_i_bar < item.second.second) { // seed_e_i_bar[e] is a descendant of the current nocde
        seed_stack.push(std::make_pair(item.first * 2 + 2, std::make_pair((item.second.first + item.second.second) / 2, item.second.second))); // Add right child
        seed_stack.push(std::make_pair(item.first * 2 + 1, std::make_pair(item.second.first, (item.second.first + item.second.second) / 2))); // Add left child
      }
      else {
        std::queue<std::pair<int, block>> seed_queue;
        seed_queue.push(std::make_pair(item.first, seed_tree[e_it][s_it++]));

        osuCrypto::PRNG prng;
        while (!seed_queue.empty()) {
          prng.SetSeed(seed_queue.front().second.b);

          if (seed_queue.front().first >= N - 1) { // Reached leaf
            partial_seeds_[e][seed_queue.front().first - (N - 1)] = seed_queue.front().second;
          }
          else { // Internal node
            seed_queue.push(std::make_pair(seed_queue.front().first * 2 + 1, prng.get<block>()));
            seed_queue.push(std::make_pair(seed_queue.front().first * 2 + 2, prng.get<block>()));
          }

          seed_queue.pop();
        }
      }
    }

    // 1.b
    osuCrypto::PRNG prng[N];
    for (auto i = 0; i < N; ++i) {
      if (i == cur_i_bar)
        continue;

      gamma_[e].resize(N);
      r_[e].resize(N);
      b_[e].resize(m);
      b_square_[e].resize(m);

      for (auto mm = 0; mm < m; ++mm){
        b_[e][mm].resize(N);
        b_square_[e][mm].resize(N);
      }

      prng[i].SetSeed(partial_seeds_[e][i].b);
      r_[e][i] = prng[i].get<block>();

      for (auto mm = 0; mm < m; ++mm) {
        b_[e][mm][i]        = NTL::ZZ_p(prng[i].get<block>().halves[0]); // NEED TO FIND A WAY TO USE ALL 128 BITS
        if (i != N - 1)
          b_square_[e][mm][i] = NTL::ZZ_p(prng[i].get<block>().halves[0]); // We ignore b_square_[e][k][N-1]
      }
    }

    // 1.c
    osuCrypto::SHA1 sha_gamma(sizeof(block));
    for (auto i = 0; i < N; ++i) {
      if (i == cur_i_bar)
        continue;

      block blk = partial_seeds_[e][i];
      sha_gamma.Reset();
      sha_gamma.Update(blk);
      if (i == N - 1) {
        unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT
        for (auto k = 0; k < m; ++k) {
          NTL::BytesFromZZ(buf, b_square[e_it][k]._ZZ_p__rep, NTL::NumBytes(b_square[e_it][k]._ZZ_p__rep));
          sha_gamma.Update(buf, NTL::NumBytes(b_square[e_it][k]._ZZ_p__rep));
        }
      }
      sha_gamma.Update(r_[e][N - 1]);
      sha_gamma.Final(gamma_[e][N - 1]);
    }

    // 1.d
    osuCrypto::SHA1 sha_h(sizeof(block));
    for (auto i = 0; i < N; ++i) {
      if (i != cur_i_bar)
        sha_h.Update(gamma_[e][i]);
      else
        sha_h.Update(gamma_i_bar[e_it]);
    }
    sha_h.Final(h_[e]); // mb need to zero it first

    // 1.e
    block gN_e, blk;
    unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT
    if (cur_i_bar != N - 1) {
      gN_e = prng[N - 1].get<block>();
      osuCrypto::SHA1 sha_omega(sizeof(block));
      for (auto i = 0; i < m; ++i) {
        NTL::BytesFromZZ(buf, s[e_it][i]._ZZ_p__rep, NTL::NumBytes(s[e_it][i]._ZZ_p__rep));
        sha_omega.Update(buf, NTL::NumBytes(s[e_it][i]._ZZ_p__rep));
      }
      sha_omega.Update(gN_e);
      sha_omega.Final(blk);
    }

    if (!eq(blk.b, omega_[e_it].b))
     return false;

    // 1.f + 1.g
    std::vector<std::vector<NTL::ZZ_p>> s_computed;
    std::vector<std::vector<NTL::ZZ_p>> alpha_computed;
    s_computed.resize(m);
    alpha_computed.resize(m);

    for (auto mm = 0; mm < m; ++mm) {
      s_computed[mm].resize(N);
      alpha_computed[mm].resize(N);
    }

    for (auto i = 0; i < N - 1; ++i) {
      if (i == cur_i_bar) {
        continue;
      }

      for (auto k = 0; k < m; ++k) {
        s_computed[k][i] = NTL::ZZ_p(prng[i].get<block>().halves[0]); // PROBABLY NEED TO CHANGE THAT
        alpha_computed[k][i] = s_computed[k][i] - b_[e][k][i];
      }
    }
    if (cur_i_bar != N - 1) {
      for (auto k = 0; k < m; ++k) {
        alpha_computed[k][N - 1] = s[e_it][k] - b_[e][k][N - 1];
      }
    }

    // 1.h
    g_[e_it] = prng_e_bar_.get<block>();
    block pi;
    osuCrypto::SHA1 sha_pi(sizeof(block));

    for (auto mm = 0; mm < m; ++mm) {
      for (auto nn = 0; nn < N; ++nn) {
        if (nn != cur_i_bar) {
          NTL::BytesFromZZ(buf, alpha_computed[mm][nn]._ZZ_p__rep, NTL::NumBytes(alpha_computed[mm][nn]._ZZ_p__rep));
          sha_pi.Update(buf, NTL::NumBytes(alpha_computed[mm][nn]._ZZ_p__rep));
        }
        else {
          NTL::BytesFromZZ(buf, alpha_i_bar[e_it][mm]._ZZ_p__rep, NTL::NumBytes(alpha_i_bar[e_it][mm]._ZZ_p__rep));
          sha_pi.Update(buf, NTL::NumBytes(alpha_i_bar[e_it][mm]._ZZ_p__rep));
        }
      }
    }
    sha_pi.Update(g_[e_it]);
    sha_pi.Final(pi);

    sha_h_pi.Update(pi); // For step 3

    // 1.i
    osuCrypto::SHA1 sha_psi(sizeof(block));

    o_[e_it].resize(N);
    for (auto i = 0; i < N; ++i) {
      if (i == cur_i_bar) {
        o_[e_it][cur_i_bar] = o_i_bar[e_it];
        continue;
      }

      for (auto l = 0; l < n; ++l) {
        NTL::ZZ_p tmp(0);

        for (auto k = 0; k < m; ++k) {
          tmp += a_[l][k] * s_computed[k][i];
        }

        o_[e_it][i] += coefficients_[e_it][l] * (t_[l] - tmp);
      }

      for (auto k = 0; k < m; ++k) {
        if (i != N - 1) {
          o_[e_it][i] += coefficients_[e_it][n + k] *
                         (alpha_computed[k][i] * (s_computed[k][i] + b_[e][k][i]) + b_square_[e][k][i] -
                          s_computed[k][i]);
        }
        else {
          o_[e_it][i] += coefficients_[e_it][n + k] *
                         (alpha_computed[k][i] * (s_computed[k][i] + b_[e][k][i]) + b_square[e_it][k] -
                          s_computed[k][i]);
        }

        NTL::BytesFromZZ(buf, o_[e_it][i]._ZZ_p__rep, NTL::NumBytes(o_[e_it][i]._ZZ_p__rep));
        sha_psi.Update(buf, NTL::NumBytes(o_[e_it][i]._ZZ_p__rep));
      }

    }

    // 1.j
    w_[e_it] = prng_e_bar_.get<block>();
    NTL::ZZ_p sigma_o(0);

    for (auto i = 0; i < N; ++i) {
      NTL::BytesFromZZ(buf, o_[e_it][i]._ZZ_p__rep, NTL::NumBytes(o_[e_it][i]._ZZ_p__rep));
      sha_psi.Update(buf, NTL::NumBytes(o_[e_it][i]._ZZ_p__rep));

      sigma_o += o_[e_it][i]; // for step 1.k
    }
    sha_psi.Update(w_[e_it]);
    sha_psi.Final(psi_[e_it]);

    sha_h_psi.Update(psi_[e_it]);

    // 1.k
    if (sigma_o != 0)
      return false;

    e_it++;
  }

  return true;

  partial_seeds = partial_seeds_;

  // 2
  osuCrypto::SHA1 sha_h_gamma(sizeof(block));
  for (auto e = 0; e < M; ++e) {
    sha_h_gamma.Update(h_[e]);
  }
  block h_gamma_computed;
  sha_h_gamma.Final(h_gamma_computed);

  if (!eq(h_gamma_.b, h_gamma_computed.b))
    return false;

  // 3
  block pi;
  sha_h_pi.Final(pi);
  if (!eq(pi.b, h_pi_.b))
    return false;

  // 4
  block psi;
  sha_h_psi.Final(psi);
  if (!eq(psi.b, h_psi_.b))
    return false;

  return true;
}