#include "prover.h"
#include "seedtree.h"

#include <cryptoTools/Crypto/sha1.h>
#include <stack>

using namespace lzkp;

Prover::Prover(const Settings &s, const NTL::Mat<NTL::ZZ_p> &a, const NTL::Vec<NTL::ZZ_p> &t, const NTL::Vec<NTL::ZZ_p> &secret)
  : M(s.M), N(s.N), q(s.q), m(s.m), n(s.n), tau(s.tau) {
  // MOVE THIS LINE TO THE DRIVER...
  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

  //a_.SetDims(n, m); // Fill with values
  //t_.SetLength(n); // Fill with values
  a_ = a;
  t_ = t;
  secret_ = secret; // Need to set private vector
}

Prover::~Prover() {
}

void Prover::r1(block &h_gamma) {
  master_seed_.resize(M);
  seed_tree_.resize(M);
  r_.resize(M);
  b_.resize(M);
  b_square_.resize(M);
  gamma_.resize(M);
  h_.resize(M);

  osuCrypto::PRNG prng;
  osuCrypto::SHA1 sha_h_gamma(sizeof(block));

  for (auto e = 0; e < M; ++e) {
    gamma_[e].resize(N);

    // 1.a
    master_seed_[e].b = osuCrypto::sysRandomSeed();
    seed_tree_[e].resize(N);
    seed_tree_[e].generate(master_seed_[e]);

    // 1.b
    r_[e].resize(N);
    b_[e].resize(m);
    b_square_[e].resize(m);

    for (auto mm = 0; mm < m; ++mm){
      b_[e][mm].resize(N);
      b_square_[e][mm].resize(N);
    }

    for (auto i = 0; i < N - 1; ++i) {
      r_[e][i] = seed_tree_[e].getBlock(i);

      for (auto mm = 0; mm < m; ++mm) {
        b_[e][mm][i]        = NTL::ZZ_p(seed_tree_[e].getBlock(i).halves[0]); // NEED TO FIND A WAY TO USE ALL 128 BITS, now working only up to 64bit
        b_square_[e][mm][i] = NTL::ZZ_p(seed_tree_[e].getBlock(i).halves[0]);
      }
    }

    // 1.c
    r_[e][N - 1] = seed_tree_[e].getBlock(N - 1);

    for (auto mm = 0; mm < m; ++mm) {
      b_[e][mm][N - 1] = NTL::ZZ_p(seed_tree_[e].getBlock(N - 1).halves[0]);
    }

    // 1.d
    for (auto k = 0; k < m; ++k) {
      b_square_[e][k][N - 1] = NTL::ZZ_p(0);
      for (auto i = 0; i < N; ++i) {
        b_square_[e][k][N - 1] += b_[e][k][i]; // Accumulate b_e,k
      }
      b_square_[e][k][N - 1] *= b_square_[e][k][N - 1]; // Calculate b_e,k^2
      for (auto i = 0; i < N - 1; ++i) {
        b_square_[e][k][N - 1] -= b_square_[e][k][i]; // Subtract b^2_e,k,i
      }
    }

    // 1.e
    osuCrypto::SHA1 sha_gamma(sizeof(block));
    for (auto i = 0; i < N; ++i) {
      block blk = seed_tree_[e].getSeed(i);
      sha_gamma.Reset(); // Calculate com(state_e,i , r_e,i) == com(seed_e,i , r_e,i)
      sha_gamma.Update(blk); // Hash seed_e,i
      if (i == N - 1) {
        unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT
        for (auto k = 0; k < m; ++k) {
          NTL::BytesFromZZ(buf, b_square_[e][k][N - 1]._ZZ_p__rep, NTL::NumBytes(b_square_[e][k][N - 1]._ZZ_p__rep));
          sha_gamma.Update(buf, NTL::NumBytes(b_square_[e][k][N - 1]._ZZ_p__rep));
        }
      }
      sha_gamma.Update(r_[e][i]); // Hash r_e,i
      sha_gamma.Final(gamma_[e][i]);
      //std::cout << "GAmma " << e << " " << i << " " << gamma_[e][i].halves[0] << " " << gamma_[e][i].halves[1] << std::endl;
    }

    // 1.f
    osuCrypto::SHA1 sha_h(sizeof(block));
    for (auto i = 0; i < N; ++i) {
      sha_h.Update(gamma_[e][i]);
    }
    sha_h.Final(h_[e]); // mb need to zero it first
    //std::cout << "H " << e << " " << h_[e].halves[0] << " " << h_[e].halves[1] << std::endl;

    sha_h_gamma.Update(h_[e]);
  }

  sha_h_gamma.Final(h_gamma_.bytes);

  h_gamma = h_gamma_;
}

void Prover::r3(const std::vector<bool> &E, std::vector<block> &seed, std::vector<block> &omegaN, block &h_pi) {
  E_ = E;

  // 1
  seed_e_bar_.b = osuCrypto::sysRandomSeed();
  prng_e_bar_.SetSeed(seed_e_bar_.b);

  // 2
  g_.resize(M);
  s_.resize(M); // We only need e_bar elements, but...
  gN_.resize(M);
  alpha_.resize(M);
  alpha_sum_.resize(M);

  osuCrypto::SHA1 sha_h_pi(sizeof(block)); // For step 3

  for (auto e = 0; e < M; ++e) {
    if (E[e]) {
      seed_.push_back(master_seed_[e]);
      continue;
    }

    // 2.a
    g_[e] = prng_e_bar_.get<block>();

    // 2.b
    s_[e].resize(m);
    alpha_[e].resize(m);
    alpha_sum_[e].resize(m);
    for (auto mm = 0; mm < m; ++mm) {
      s_[e][mm].resize(N);
      alpha_[e][mm].resize(N);
    }

    for (auto i = 0; i < N - 1; ++i) {
      for (auto mm = 0; mm < m; ++mm) {
        s_[e][mm][i] = NTL::ZZ_p(seed_tree_[e].getBlock(i).halves[0]); // PROBABLY NEED TO CHANGE THAT
      }
    }

    // 2.c
    gN_[e] = seed_tree_[e].getBlock(N - 1);

    // 2.d
    for (auto k = 0; k < m; ++k) {
      s_[e][k][N - 1] = secret_[k];

      for (auto i = 0; i < N - 1; ++i) {
        s_[e][k][N - 1] -= s_[e][k][i];
      }
    }

    // 2.e + 2.h
    for (auto k = 0; k < m; ++k) {
      alpha_sum_[e][k] = NTL::ZZ_p(0);

      for (auto i = 0; i < N; ++i) {
        alpha_[e][k][i] = s_[e][k][i] - b_[e][k][i];
        alpha_sum_[e][k] += alpha_[e][k][i];
        //std::cout << "alpha " << k << " " << i << " " << alpha_[e][k][i] << std::endl;
      }
      //std::cout << "alpha_sum " << e << " " << k << " " << alpha_sum_[e][k] << std::endl;
    }

    // 2.f
    block blk;
    osuCrypto::SHA1 sha_omegaN(sizeof(block));
    unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT
    for (auto i = 0; i < m; ++i) {
      NTL::BytesFromZZ(buf, s_[e][i][N - 1]._ZZ_p__rep, NTL::NumBytes(s_[e][i][N - 1]._ZZ_p__rep));
      sha_omegaN.Update(buf, NTL::NumBytes(s_[e][i][N - 1]._ZZ_p__rep));
    }
    sha_omegaN.Update(gN_[e]);
    sha_omegaN.Final(blk);

    omegaN_.push_back(blk);

    // 2.g
    block pi;
    osuCrypto::SHA1 sha_pi(sizeof(block));

    for (auto mm = 0; mm < m; ++mm) {
      for (auto nn = 0; nn < N; ++nn) {
        NTL::BytesFromZZ(buf, alpha_[e][mm][nn]._ZZ_p__rep, NTL::NumBytes(alpha_[e][mm][nn]._ZZ_p__rep));
        sha_pi.Update(buf, NTL::NumBytes(alpha_[e][mm][nn]._ZZ_p__rep));
      }
    }
    sha_pi.Update(g_[e]);
    sha_pi.Final(pi);

    //std::cout << "PI E " << e << " " << pi.halves[0] << pi.halves[1] << std::endl;

    sha_h_pi.Update(pi);
  }

  sha_h_pi.Final(h_pi_);

  seed = seed_;
  omegaN = omegaN_;
  h_pi = h_pi_;
}

void Prover::r5(const block &seed_ell, block &h_psi) {
  seed_ell_ = seed_ell;

  prng_seed_ell_.SetSeed(seed_ell_.b);
  coefficients_.resize(M - tau);

  o_.resize(M - tau);
  psi_.resize(M - tau);
  w_.resize(M - tau);

  int e_it = 0;
  osuCrypto::SHA1 sha_h_psi(sizeof(block));

  for (auto e = 0; e < M; ++e) {
    if (E_[e]) {
      continue;
    }

    osuCrypto::SHA1 sha_psi(sizeof(block));
    unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT

    // 1
    coefficients_[e_it].resize(n + m);

    for (auto i = 0; i < n; ++i) {
      coefficients_[e_it][i] = NTL::ZZ_p(prng_seed_ell_.get<block>().halves[0]);
    }
    for (auto i = 0; i < m; ++i) {
      coefficients_[e_it][n + i] = NTL::ZZ_p(prng_seed_ell_.get<block>().halves[0]);
    }

    // 2.a
    o_[e_it].resize(N);
    for (auto i = 0; i < N; ++i) {
      o_[e_it][i] = NTL::ZZ_p(0);

      for (auto l = 0; l < n; ++l) {
        NTL::ZZ_p tmp(0);

        for (auto k = 0; k < m; ++k) {
          tmp += a_[l][k] * s_[e][k][i];
        }

        o_[e_it][i] += coefficients_[e_it][l] * ((t_[l] / N) - tmp);
      }

      for (auto k = 0; k < m; ++k) {
        //o_[e_it][i] += coefficients_[e_it][n + k] *
//                       (alpha_sum_[e][k] * (s_[e][k][i] + b_[e][k][i]) + b_square_[e][k][i] - s_[e][k][i]);
        o_[e_it][i] += coefficients_[e_it][n + k] *
                       (alpha_sum_[e][k] * (s_[e][k][i] + b_[e][k][i]) + b_square_[e][k][i] - s_[e][k][i]);
      }

      NTL::BytesFromZZ(buf, o_[e_it][i]._ZZ_p__rep, NTL::NumBytes(o_[e_it][i]._ZZ_p__rep));
      sha_psi.Update(buf, NTL::NumBytes(o_[e_it][i]._ZZ_p__rep)); // For step 2.b

      //std::cout << "O " << e << " " << i << " " << o_[e_it][i] << std::endl;
    }

    // 2.b
    w_[e_it] = prng_e_bar_.get<block>();
    sha_psi.Update(w_[e_it]);
    sha_psi.Final(psi_[e_it]);

    sha_h_psi.Update(psi_[e_it]);

    //std::cout << "WE " << e << " " << w_[e_it].halves[0] << " " << w_[e_it].halves[1] << std::endl;
    //std::cout << "PSI " << e << " " << psi_[e_it].halves[0] << psi_[e_it].halves[1] << std::endl;

    e_it++;
  }

  sha_h_psi.Final(h_psi_);

  h_psi = h_psi_;
}

void Prover::r7(const std::vector<int> &i_bar, block &seed_e_bar, std::vector<std::vector<block>> &seed_tree,
                std::vector<block> &gamma_i_bar, std::vector<std::vector<NTL::ZZ_p>> &alpha_i_bar, std::vector<NTL::ZZ_p> &o_i_bar,
                std::vector<std::vector<NTL::ZZ_p>> &b_square, std::vector<std::vector<NTL::ZZ_p>> &s) {
  seed_e_bar = seed_e_bar_;
  seed_tree.resize(M - tau);
  gamma_i_bar.resize(M - tau);
  alpha_i_bar.resize(M - tau);
  o_i_bar.resize(M - tau);
  b_square.resize(M - tau);
  s.resize(M - tau);

  int e_it = 0;

  for (auto e = 0; e < M; ++e) {
    if (E_[e]) {
      continue;
    }

    const auto cur_i_bar = i_bar[e_it];

    seed_tree[e_it].clear();

    std::stack<std::pair<int, std::pair<int, int>>> seed_stack;

    seed_stack.push(std::make_pair(0, std::make_pair(0, N)));

    while (!seed_stack.empty()) {
      auto item = seed_stack.top();
      seed_stack.pop();

      // Check for leaf
      if (item.second.second - item.second.first == 1) {
        if (item.second.first != cur_i_bar) {
          seed_tree[e_it].push_back(seed_tree_[e].getSeed(item.second.first));
        }
        continue;
      }

      if (cur_i_bar >= item.second.first && cur_i_bar < item.second.second) { // seed_e_i_bar[e] is a descendant of the current nocde
        seed_stack.push(std::make_pair(item.first * 2 + 2, std::make_pair((item.second.first + item.second.second) / 2, item.second.second))); // Add right child
        seed_stack.push(std::make_pair(item.first * 2 + 1, std::make_pair(item.second.first, (item.second.first + item.second.second) / 2))); // Add left child
      }
      else {
        seed_tree[e_it].push_back(seed_tree_[e][item.first]);
      }
    }

    gamma_i_bar[e_it] = gamma_[e][cur_i_bar];

    alpha_i_bar[e_it].resize(m);
    for (auto k = 0; k < m; ++k) {
      alpha_i_bar[e_it][k] = alpha_[e][k][cur_i_bar];
    }

    o_i_bar[e_it] = o_[e_it][cur_i_bar];

    if (cur_i_bar != N - 1) {
      b_square[e_it].resize(m);
      s[e_it].resize(m);

      for (auto k = 0; k < m; ++k) {
        b_square[e_it][k] = b_square_[e][k][N - 1];
        s[e_it][k] = s_[e][k][N - 1];
      }
    }

    e_it++;
  }
}