#ifndef __LZKP_SAC_PROVER_H_FILE___
#define __LZKP_SAC_PROVER_H_FILE___


#include <stack>
#include <cryptoTools/Crypto/sha1.h>

#include "seedtree.h"
#include "parameters.h"



namespace lzkp {


template <class FieldType>
class SacProver {
public:
  SacProver(const Parameters &s, const std::vector<std::vector<FieldType>> &a, const std::vector<FieldType> &t, const std::vector<FieldType> &secret);
  ~SacProver();

  void r1(const block &master_seed);
  void r3(); // Local variables seed_ell_, seed_global_, g_, w_, u_ must be set before calling this method
  void r5(std::vector<block> &seed_tree, block &gamma_i_bar, std::vector<FieldType> &alpha_i_bar,
          FieldType &o_i_bar, FieldType &v_i_bar, std::vector<FieldType> &b_square, std::vector<FieldType> &s, std::vector<FieldType> &s_square);
             // Local variable i_bar_ must be set before calling this method

private:
public:
  // Public known values
  const std::vector<std::vector<FieldType>> &a_;
  const std::vector<FieldType> &t_;

  // CacProver's secret
  const std::vector<FieldType> secret_;

  const int N;
  const int n;
  const int m;

  // Local values
  block master_seed_;
  SeedTree seed_tree_;
  std::vector<block> r_;
  std::vector<std::vector<FieldType>> s_, s_square_;
  std::vector<std::vector<FieldType>> b_, b_square_;
  std::vector<block> gamma_;
  block h_;

  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;
  std::vector<FieldType> ep_, be_, ga_;
  block g_;
  std::vector<std::vector<FieldType>> alpha_;
  std::vector<FieldType> alpha_sum_;
  block pi_;
  std::vector<FieldType> o_;
  block w_;
  block psi_;
  block seed_global_;
  osuCrypto::PRNG prng_seed_global_;
  std::vector<FieldType> de_;
  std::vector<FieldType> v_;
  block u_;
  block theta_;

  int i_bar_;

//  block gN_; // g_{e,N}
//  block omegaN_;

//
//  std::vector<FieldType> coefficients_;

//  block w_;
//
};

template <class FieldType>
SacProver<FieldType>::SacProver(const Parameters &s, const std::vector<std::vector<FieldType>> &a, const std::vector<FieldType> &t, const std::vector<FieldType> &secret)
    : a_(a), t_(t), secret_(secret), N(s.N), n(s.n), m(s.m) {
}

template <class FieldType>
SacProver<FieldType>::~SacProver() {
}

template <class FieldType>
void SacProver<FieldType>::r1(const block &master_seed) {
  master_seed_ = master_seed;

  // 1.a
  seed_tree_.resize(N);
  seed_tree_.generate(master_seed_);

  // 1.b
  r_.resize(N);
  s_.resize(m);
  s_square_.resize(m);
  b_.resize(m);
  b_square_.resize(m);

  for (auto mm = 0; mm < m; ++mm) {
    s_[mm].resize(N);
    s_square_[mm].resize(N);
    b_[mm].resize(N);
    b_square_[mm].resize(N);
  }

  for (auto i = 0; i < N - 1; ++i) {
    r_[i] = seed_tree_.getBlock(i);

    for (auto mm = 0; mm < m; ++mm) {
      s_[mm][i]        = FieldType(seed_tree_.getBlock(i).halves[0]); // NEED TO FIND A WAY TO USE ALL 128 BITS, now working only up to 64bit
      s_square_[mm][i] = FieldType(seed_tree_.getBlock(i).halves[0]);
      b_[mm][i]        = FieldType(seed_tree_.getBlock(i).halves[0]); // NEED TO FIND A WAY TO USE ALL 128 BITS, now working only up to 64bit
      b_square_[mm][i] = FieldType(seed_tree_.getBlock(i).halves[0]);
    }
  }

  // 1.c
  r_[N - 1] = seed_tree_.getBlock(N - 1);

  for (auto mm = 0; mm < m; ++mm) {
    b_[mm][N - 1] = FieldType(seed_tree_.getBlock(N - 1).halves[0]);
  }

  // 1.d
  for (auto k = 0; k < m; ++k) {
    b_square_[k][N - 1] = FieldType(0);
    for (auto i = 0; i < N; ++i) {
      b_square_[k][N - 1] += b_[k][i]; // Accumulate b_k
    }
    b_square_[k][N - 1] *= b_square_[k][N - 1]; // Calculate (b_k)^2
    for (auto i = 0; i < N - 1; ++i) {
      b_square_[k][N - 1] -= b_square_[k][i]; // Subtract b^2_e,k,i
    }
  }

  // 1.e - Check if this step is needed at all
  for (auto k = 0; k < m; ++k) {
    s_[k][N - 1] = secret_[k];
    s_square_[k][N - 1] = secret_[k] * secret_[k];

    for (auto i = 0; i < N - 1; ++i) {
      s_[k][N - 1] -= s_[k][i];
      s_square_[k][N - 1] -= s_square_[k][i];
    }
  }

  // 1.f
  gamma_.resize(N);
  osuCrypto::SHA1 sha_gamma(sizeof(block));

  for (auto i = 0; i < N; ++i) {
    block blk = seed_tree_.getSeed(i);
    sha_gamma.Reset(); // Calculate com(state_e,i , r_e,i) == com(seed_e,i , r_e,i)
    sha_gamma.Update(blk); // Hash seed_e,i
    if (i == N - 1) {
      for (auto k = 0; k < m; ++k) { // TODO: optimize
        sha_gamma.Update(s_[k][N - 1].elem); // TODO: check POS
        sha_gamma.Update(s_square_[k][N - 1].elem);
        sha_gamma.Update(b_square_[k][N - 1].elem);
      }
    }
    sha_gamma.Update(r_[i]); // Hash r_e,i
    sha_gamma.Final(gamma_[i]);
  }

  // 1.g
  osuCrypto::SHA1 sha_h(sizeof(block));
  for (auto i = 0; i < N; ++i) {
    sha_h.Update(gamma_[i]);
  }
  sha_h.Final(h_); // mb need to zero it first
}

template <class FieldType>
void SacProver<FieldType>::r3() {
  // 1
  ep_.resize(m);
  be_.resize(n);
  ga_.resize(m);

  prng_seed_ell_.SetSeed(seed_ell_.b);

  for (auto i = 0; i < m; ++i) {
    ep_[i] = FieldType(prng_seed_ell_.get<block>().halves[0]);
    ga_[i] = FieldType(prng_seed_ell_.get<block>().halves[0]);
  }

  for (auto i = 0; i < n; ++i) {
    be_[i] = FieldType(prng_seed_ell_.get<block>().halves[0]);
  }

  // 3.a + 3.b + 3.c (combined in loop)
  alpha_.resize(m);
  alpha_sum_.resize(m);

  for (auto mm = 0; mm < m; ++mm) {
    alpha_[mm].resize(N);
    alpha_sum_[mm] = FieldType(0); // For step 3.b
  }

  osuCrypto::SHA1 sha_pi(sizeof(block)); // For step 3.c

  for (auto i = 0; i < N; ++i) {
    for (auto k = 0; k < m; ++k) {
      alpha_[k][i] = s_[k][i] - ep_[k] * b_[k][i];
      alpha_sum_[k] += alpha_[k][i]; // 3.b

      sha_pi.Update(alpha_[k][i].elem);
    }
  }
  sha_pi.Update(g_);
  sha_pi.Final(pi_);

  // 3.d + 3.e
  o_.resize(N);

  osuCrypto::SHA1 sha_psi(sizeof(block));

  for (auto i = 0; i < N; ++i) {
    o_[i] = FieldType(0);

    for (auto l = 0; l < n; ++l) {
      FieldType tmp(0);

      for (auto k = 0; k < m; ++k) {
        tmp += a_[l][k] * s_[k][i];
      }

      o_[i] += be_[l] * ((t_[l] / FieldType(N)) - tmp);
    }

    for (auto k = 0; k < m; ++k) {
      o_[i] += ga_[k] * (s_square_[k][i] - s_[k][i]);
    }

    sha_psi.Update(o_[i].elem); // For step 2.b
  }

  sha_psi.Update(w_);
  sha_psi.Final(psi_);

  // 3.f + 3.g
  v_.resize(N);

  de_.resize(m);
  prng_seed_global_.SetSeed(seed_global_.b);

  for (auto k = 0; k < m; ++k) {
    de_[k] = FieldType(prng_seed_global_.get<block>().halves[0]);
  }

  osuCrypto::SHA1 sha_theta(sizeof(block));

  for (auto i = 0; i < N; ++i) {
    v_[i] = FieldType(0);

    for (auto k = 0; k < m; ++k) {
      v_[i] += de_[k] * (s_square_[k][i] - alpha_sum_[k] * (s_[k][i] + ep_[k] * b_[k][i]) - ((ep_[k] * ep_[k]) * b_square_[k][i]));
    }

    sha_theta.Update(v_[i].elem); // For step 2.b
  }

  sha_theta.Update(u_);
  sha_theta.Final(theta_);
}

template <class FieldType>
void SacProver<FieldType>::r5(std::vector<block> &seed_tree, block &gamma_i_bar, std::vector<FieldType> &alpha_i_bar,
              FieldType &o_i_bar, FieldType &v_i_bar, std::vector<FieldType> &b_square, std::vector<FieldType> &s, std::vector<FieldType> &s_square) {
  // 1
  seed_tree.clear();

  std::stack<std::pair<int, std::pair<int, int>>> seed_stack;

  seed_stack.push(std::make_pair(0, std::make_pair(0, N)));

  while (!seed_stack.empty()) {
    auto item = seed_stack.top();
    seed_stack.pop();

    // Check for leaf
    if (item.second.second - item.second.first == 1) {
      if (item.second.first != i_bar_) {
        seed_tree.push_back(seed_tree_.getSeed(item.second.first));
      }
      continue;
    }

    if (i_bar_ >= item.second.first && i_bar_ < item.second.second) { // seed_e_i_bar[e] is a descendant of the current node
      seed_stack.push(std::make_pair(item.first * 2 + 2, std::make_pair((item.second.first + item.second.second) / 2, item.second.second))); // Add right child
      seed_stack.push(std::make_pair(item.first * 2 + 1, std::make_pair(item.second.first, (item.second.first + item.second.second) / 2))); // Add left child
    }
    else {
      seed_tree.push_back(seed_tree_[item.first]);
    }
  }

  gamma_i_bar = gamma_[i_bar_];

  alpha_i_bar.resize(m);
  for (auto k = 0; k < m; ++k) {
    alpha_i_bar[k] = alpha_[k][i_bar_];
  }

  o_i_bar = o_[i_bar_];

  v_i_bar = v_[i_bar_];

  if (i_bar_ != N - 1) {
    b_square.resize(m);
    s.resize(m);
    s_square.resize(m);

    for (auto k = 0; k < m; ++k) {
      b_square[k] = b_square_[k][N - 1];
      s[k] = s_[k][N - 1];
      s_square[k] = s_square_[k][N - 1];
    }
  }
}


}


#endif // __LZKP_SAC_PROVER_H_FILE___