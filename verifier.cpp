//
// Created by roee on 7/26/18.
//

#include "verifier.h"

#include <stack>
#include <queue>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Crypto/sha1.h>

#include "seedtree.h"


using namespace lzkp;

Verifier::Verifier(const Settings &s, const NTL::Mat<NTL::ZZ_p> &a, const NTL::Vec<NTL::ZZ_p> &t)
    : N(s.N), q(s.q), m(s.m), n(s.n), a_(a), t_(t) {
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***
}

Verifier::~Verifier() {
}

void Verifier::r4() {
//  coefficients_.resize(M - tau);

//  seed_tree_.resize(M);
//  r_.resize(M);
//  b_.resize(M);
//  b_square_.resize(M);
//  gamma_.resize(M);
//  h_.resize(M);

  seed_tree_.resize(N);
  seed_tree_.generate(seed_); // Generate seed_tree
  gamma_.resize(N);

  r_.resize(N);
  b_.resize(m);
  b_square_.resize(m);

  for (auto mm = 0; mm < m; ++mm){
    b_[mm].resize(N);
    b_square_[mm].resize(N);
  }

  // *1* - 1.b
  for (auto i = 0; i < N - 1; ++i) {
    r_[i] = seed_tree_.getBlock(i);

    for (auto mm = 0; mm < m; ++mm) {
      b_[mm][i]        = NTL::ZZ_p(seed_tree_.getBlock(i).halves[0]); // NEED TO FIND A WAY TO USE ALL 128 BITS
      b_square_[mm][i] = NTL::ZZ_p(seed_tree_.getBlock(i).halves[0]);
    }
  }

  // *1* - 1.c
  r_[N - 1] = seed_tree_.getBlock(N - 1);

  for (auto mm = 0; mm < m; ++mm) {
    b_[mm][N - 1] = NTL::ZZ_p(seed_tree_.getBlock(N - 1).halves[0]);
  }

  // *1* - 1.d
  for (auto k = 0; k < m; ++k) {
    b_square_[k][N - 1] = NTL::ZZ_p(0);
    for (auto i = 0; i < N; ++i) {
      b_square_[k][N - 1] += b_[k][i];
    }
    b_square_[k][N - 1] *= b_square_[k][N - 1];
    for (auto i = 0; i < N - 1; ++i) {
      b_square_[k][N - 1] -= b_square_[k][i];
    }
  }

  // *1* - 1.e
  osuCrypto::SHA1 sha_gamma(sizeof(block));
  for (auto i = 0; i < N - 1; ++i) {
    block blk = seed_tree_.getSeed(i);
    sha_gamma.Reset();
    sha_gamma.Update(blk);
    sha_gamma.Update(r_[i]);
    sha_gamma.Final(gamma_[i]);
  }
  block blk = seed_tree_.getSeed(N - 1);
  sha_gamma.Reset();
  sha_gamma.Update(blk);
  unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT
  for (auto i = 0; i < m; ++i) {
    NTL::BytesFromZZ(buf, b_square_[i][N - 1]._ZZ_p__rep, NTL::NumBytes(b_square_[i][N - 1]._ZZ_p__rep));
    sha_gamma.Update(buf, NTL::NumBytes(b_square_[i][N - 1]._ZZ_p__rep));
  }
  sha_gamma.Update(r_[N - 1]);
  sha_gamma.Final(gamma_[N - 1]);

  // *1* - 1.f
  osuCrypto::SHA1 sha_h(sizeof(block));
  for (auto i = 0; i < N; ++i) {
    sha_h.Update(gamma_[i]);
  }
  sha_h.Final(h_); // mb need to zero it first
  //std::cout << "Vr4 H " << e << " " << h_[e].halves[0] << " " << h_[e].halves[1] << std::endl;
}

bool Verifier::r8(const std::vector<block> &seed_tree, const block &gamma_i_bar,
                  const std::vector<NTL::ZZ_p> &alpha_i_bar, const NTL::ZZ_p &o_i_bar,
                  const std::vector<NTL::ZZ_p> &b_square, const std::vector<NTL::ZZ_p> &s) {
  reject_ = false;

  // 1.a
  partial_seeds_.resize(N);
  int s_it = 0;

  std::stack<std::pair<int, std::pair<int, int>>> seed_stack;

  seed_stack.push(std::make_pair(0, std::make_pair(0, N)));

  while (!seed_stack.empty()) {
    auto item = seed_stack.top();
    seed_stack.pop();

    // Check for leaf
    if (item.second.second - item.second.first == 1) {
      if (item.second.first != i_bar_) {
        partial_seeds_[item.second.first] = seed_tree[s_it++];
      }
      continue;
    }

    if (i_bar_ >= item.second.first &&
        i_bar_ < item.second.second) { // seed_e_i_bar[e] is a descendant of the current nocde
      seed_stack.push(std::make_pair(item.first * 2 + 2, std::make_pair((item.second.first + item.second.second) / 2,
                                                                        item.second.second))); // Add right child
      seed_stack.push(std::make_pair(item.first * 2 + 1, std::make_pair(item.second.first,
                                                                        (item.second.first + item.second.second) /
                                                                        2))); // Add left child
    } else {
      std::queue<std::pair<int, block>> seed_queue;
      seed_queue.push(std::make_pair(item.first, seed_tree[s_it++]));

      osuCrypto::PRNG prng;
      while (!seed_queue.empty()) {
        prng.SetSeed(seed_queue.front().second.b);

        if (seed_queue.front().first >= N - 1) { // Reached leaf
          partial_seeds_[seed_queue.front().first - (N - 1)] = seed_queue.front().second;
        } else { // Internal node
          seed_queue.push(std::make_pair(seed_queue.front().first * 2 + 1, prng.get<block>()));
          seed_queue.push(std::make_pair(seed_queue.front().first * 2 + 2, prng.get<block>()));
        }

        seed_queue.pop();
      }
    }
  }

  // 1.b
  std::vector<osuCrypto::PRNG> prng(N);
  for (auto i = 0; i < N; ++i) {
    if (i == i_bar_)
      continue;

    gamma_.resize(N);
    r_.resize(N);
    b_.resize(m);
    b_square_.resize(m);

    for (auto mm = 0; mm < m; ++mm){
      b_[mm].resize(N);
      b_square_[mm].resize(N);
    }

    prng[i].SetSeed(partial_seeds_[i].b);
    r_[i] = prng[i].get<block>();

    for (auto mm = 0; mm < m; ++mm) {
      b_[mm][i]        = NTL::ZZ_p(prng[i].get<block>().halves[0]); // NEED TO FIND A WAY TO USE ALL 128 BITS
      if (i != N - 1)
        b_square_[mm][i] = NTL::ZZ_p(prng[i].get<block>().halves[0]); // We ignore b_square_[e][k][N-1]
    }
  }

  // 1.c
  osuCrypto::SHA1 sha_gamma(sizeof(block));
  for (auto i = 0; i < N; ++i) {
    if (i == i_bar_)
      continue;

    block blk = partial_seeds_[i];
    sha_gamma.Reset();
    sha_gamma.Update(blk);
    if (i == N - 1) {
      unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT
      for (auto k = 0; k < m; ++k) {
        NTL::BytesFromZZ(buf, b_square[k]._ZZ_p__rep, NTL::NumBytes(b_square[k]._ZZ_p__rep));
        sha_gamma.Update(buf, NTL::NumBytes(b_square[k]._ZZ_p__rep));
      }
    }
    sha_gamma.Update(r_[i]);
    sha_gamma.Final(gamma_[i]);

  }

  // 1.d
  osuCrypto::SHA1 sha_h(sizeof(block));
  for (auto i = 0; i < N; ++i) {
    if (i != i_bar_)
      sha_h.Update(gamma_[i]);
    else
      sha_h.Update(gamma_i_bar);
  }
  sha_h.Final(h_); // mb need to zero it first
  //std::cout << "Vr8 H " << e << " " << h_[e].halves[0] << " " << h_[e].halves[1] << std::endl;

  // 1.e
  unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT
  if (i_bar_ != N -1) {
    block gN_e, blk;
    if (i_bar_ != N - 1) {
      gN_e = prng[N - 1].get<block>();
      osuCrypto::SHA1 sha_omega(sizeof(block));
      for (auto i = 0; i < m; ++i) {
        NTL::BytesFromZZ(buf, s[i]._ZZ_p__rep, NTL::NumBytes(s[i]._ZZ_p__rep));
        sha_omega.Update(buf, NTL::NumBytes(s[i]._ZZ_p__rep));
      }
      sha_omega.Update(gN_e);
      sha_omega.Final(blk);
    }

    if (!eq(blk.b, omegaN_.b)) {
//      std::cout << "1\t" << blk.b << " " << omegaN_.b << std::endl;

      reject_ = true;
      return false;
    }
  }

  // 1.f + 1.g + 1.i
  std::vector<std::vector<NTL::ZZ_p>> s_computed;
  std::vector<std::vector<NTL::ZZ_p>> alpha_computed;
  std::vector<NTL::ZZ_p> alpha_sum_computed;
  s_computed.resize(m);
  alpha_computed.resize(m);
  alpha_sum_computed.resize(m);

  for (auto mm = 0; mm < m; ++mm) {
    s_computed[mm].resize(N);
    alpha_computed[mm].resize(N);
    alpha_sum_computed[mm] = NTL::ZZ_p(0);
  }

  // 1.f
  for (auto i = 0; i < N - 1; ++i) {
    if (i == i_bar_)
      continue;

    for (auto k = 0; k < m; ++k)
      s_computed[k][i] = NTL::ZZ_p(prng[i].get<block>().halves[0]); // PROBABLY NEED TO CHANGE THAT
  }

  // 1.g + 1.i
  for (auto i = 0; i < N; ++i) {
    for (auto k = 0; k < m; ++k) {
      if (i != i_bar_) {
        if (i != N - 1)
          alpha_computed[k][i] = s_computed[k][i] - b_[k][i];
        else
          alpha_computed[k][i] = s[k] - b_[k][i];

        alpha_sum_computed[k] += alpha_computed[k][i]; // 1.i
      }
      else {
        alpha_sum_computed[k] += alpha_i_bar[k];
      }
    }
  }

  // 1.h
  osuCrypto::SHA1 sha_pi(sizeof(block));

  for (auto mm = 0; mm < m; ++mm) {
    for (auto nn = 0; nn < N; ++nn) {
      if (nn != i_bar_) {
        NTL::BytesFromZZ(buf, alpha_computed[mm][nn]._ZZ_p__rep, NTL::NumBytes(alpha_computed[mm][nn]._ZZ_p__rep));
        sha_pi.Update(buf, NTL::NumBytes(alpha_computed[mm][nn]._ZZ_p__rep));
      }
      else {
        NTL::BytesFromZZ(buf, alpha_i_bar[mm]._ZZ_p__rep, NTL::NumBytes(alpha_i_bar[mm]._ZZ_p__rep));
        sha_pi.Update(buf, NTL::NumBytes(alpha_i_bar[mm]._ZZ_p__rep));
      }
    }
  }
  sha_pi.Update(g_);
  sha_pi.Final(pi_);

  // 1.j
  osuCrypto::SHA1 sha_psi(sizeof(block));

  o_.resize(N);
  for (auto i = 0; i < N; ++i) {
    o_[i] = NTL::ZZ_p(0);

    if (i == i_bar_) {
      o_[i] = o_i_bar;
      //std::cout << "V O " << e << " " << i << " " << o_[e_it][i] << std::endl;

      continue;
    }

    for (auto l = 0; l < n; ++l) {
      NTL::ZZ_p tmp(0);

      for (auto k = 0; k < m; ++k) {
//        std::cout <<  s_computed[k][i] << " " << a_[l][k] <<  " " << s[k] << std::endl;
        if (i != N - 1)
          tmp += a_[l][k] * s_computed[k][i];
        else
          tmp += a_[l][k] * s[k];
      }

      o_[i] += coefficients_[l] * ((t_[l] / N) - tmp);
    }

    for (auto k = 0; k < m; ++k) {
      if (i != N - 1) {
        o_[i] += coefficients_[n + k] *
                       (alpha_sum_computed[k] * (s_computed[k][i] + b_[k][i]) + b_square_[k][i] -
                        s_computed[k][i]);
      }
      else {
        o_[i] += coefficients_[n + k] *
                       (alpha_sum_computed[k] * (s[k] + b_[k][i]) + b_square[k] -
                        s[k]);
      }
    }
    //std::cout << "V O " << e << " " << i << " " << o_[e_it][i] << std::endl;
  }

  // 1.k
  sigma_o = 0;

  for (auto i = 0; i < N; ++i) {
    NTL::BytesFromZZ(buf, o_[i]._ZZ_p__rep, NTL::NumBytes(o_[i]._ZZ_p__rep));
    sha_psi.Update(buf, NTL::NumBytes(o_[i]._ZZ_p__rep));

    sigma_o += o_[i]; // for step 1.k
  }
  sha_psi.Update(w_);
  sha_psi.Final(psi_);

  //std::cout << "V WE " << e << " " << w_[e_it].halves[0] << " " << w_[e_it].halves[1] << std::endl;
  //std::cout << "V PSI " << e << " " << psi_[e_it].halves[0] << psi_[e_it].halves[1] << std::endl;

  // 1.l
  if (sigma_o != 0) {
    //std::cout << "BAR: " << cur_i_bar << " 2\t" << sigma_o << std::endl;
    reject_ = true;
    return false;
  }


//    for (auto k = 0; k < m; ++k) {
      //std::cout << "V alpha_sum " << k << " " << alpha_sum_computed[k] << std::endl;
//    }



    //std::cout << "V PI E " << e << " " << pi.halves[0] << pi.halves[1] << std::endl;

  return true;
}
