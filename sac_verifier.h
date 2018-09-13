#ifndef __LZKP_SAC_VERIFIER_H_FILE___
#define __LZKP_SAC_VERIFIER_H_FILE___


#include <stack>
#include <queue>
#include <cryptoTools/Crypto/sha1.h>
#include <cryptoTools/Crypto/PRNG.h>

#include "seedtree.h"
#include "parameters.h"


namespace lzkp {


template <class FieldType>
class SacVerifier {
public:
  SacVerifier(const Parameters &s, FieldType **&a, FieldType *&t);
  ~SacVerifier();

  void r2(); // Local variable seed_ell_ must be set before calling this method
  bool r6(const std::vector<block> &seed_tree, const block &gamma_i_bar,
          const std::vector<FieldType> &alpha_i_bar, const FieldType &o_i_bar, const FieldType &v_i_bar,
          const std::vector<FieldType> &b_square, const std::vector<FieldType> &s, const std::vector<FieldType> &s_square);
    // Local variables seed_global_, g_, w_, u_  must be set before calling this method

private:
public:
  // Public known values
  FieldType **&a_;
  FieldType *&t_;

  const int N;
  const int n;
  const int m;

  // Local values
  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;
  std::vector<FieldType> ep_, be_, ga_;

  int i_bar_;

  SeedTree seed_tree_;
  std::vector<block> r_;
  std::vector<std::vector<FieldType>> s_, s_square_;
  std::vector<std::vector<FieldType>> b_, b_square_;
  std::vector<block> gamma_;
  block h_;

  bool reject_;
  std::vector<block> partial_seeds_;
  block g_;
  block pi_;
  std::vector<FieldType> o_;
  FieldType sigma_o;
  block w_;
  block psi_;
  block seed_global_;
  osuCrypto::PRNG prng_seed_global_;
  std::vector<FieldType> de_;
  std::vector<FieldType> v_;
  FieldType sigma_v;
  block u_;
  block theta_;
};

template <class FieldType>
SacVerifier<FieldType>::SacVerifier(const Parameters &s, FieldType **&a, FieldType *&t)
    : a_(a), t_(t), N(s.N), n(s.n), m(s.m) {
}

template <class FieldType>
SacVerifier<FieldType>::~SacVerifier() {
}

template <class FieldType>
void SacVerifier<FieldType>::r2() {
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
}

template <class FieldType>
bool SacVerifier<FieldType>::r6(const std::vector<block> &seed_tree, const block &gamma_i_bar,
                        const std::vector<FieldType> &alpha_i_bar, const FieldType &o_i_bar, const FieldType &v_i_bar,
                        const std::vector<FieldType> &b_square, const std::vector<FieldType> &s, const std::vector<FieldType> &s_square) {
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

    if (i_bar_ >= item.second.first && i_bar_ < item.second.second) { // seed_e_i_bar[e] is a descendant of the current node
      seed_stack.push(std::make_pair(item.first * 2 + 2, std::make_pair((item.second.first + item.second.second) / 2,
                                                                        item.second.second))); // Add right child
      seed_stack.push(std::make_pair(item.first * 2 + 1, std::make_pair(item.second.first,
                                                                        (item.second.first + item.second.second) /
                                                                        2))); // Add left child
    }
    else {
      std::queue<std::pair<int, block>> seed_queue;
      seed_queue.push(std::make_pair(item.first, seed_tree[s_it++]));

      osuCrypto::PRNG prng;
      while (!seed_queue.empty()) {
        prng.SetSeed(seed_queue.front().second.b);

        if (seed_queue.front().first >= N - 1) { // Reached leaf
          partial_seeds_[seed_queue.front().first - (N - 1)] = seed_queue.front().second;
        } else { // Internal node
          block bb;
          bb = prng.get<block>();
          seed_queue.push(std::make_pair(seed_queue.front().first * 2 + 1, bb));
          bb = prng.get<block>();
          seed_queue.push(std::make_pair(seed_queue.front().first * 2 + 2, bb));
        }

        seed_queue.pop();
      }
    }
  }

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

  std::vector<osuCrypto::PRNG> prng(N);
  for (auto i = 0; i < N; ++i) {
    if (i == i_bar_)
      continue;

    prng[i].SetSeed(partial_seeds_[i].b);
    r_[i] = prng[i].get<block>();

    for (auto mm = 0; mm < m; ++mm) {
      if (i != N - 1) {
        s_[mm][i] = FieldType(prng[i].get<block>().halves[0]);
        s_square_[mm][i] = FieldType(prng[i].get<block>().halves[0]);
      }
      b_[mm][i] = FieldType(prng[i].get<block>().halves[0]); // NEED TO FIND A WAY TO USE ALL 128 BITS
      if (i != N - 1)
        b_square_[mm][i] = FieldType(prng[i].get<block>().halves[0]); // We ignore b_square_[e][k][N-1]
    }
  }

  // 1.c
  gamma_.resize(N);

  osuCrypto::SHA1 sha_gamma(sizeof(block));
  for (auto i = 0; i < N; ++i) {
    if (i == i_bar_)
      continue;

    block blk = partial_seeds_[i];
    sha_gamma.Reset();
    sha_gamma.Update(blk);
    if (i == N - 1) {
      for (auto k = 0; k < m; ++k) {
        sha_gamma.Update(s[k].elem); // TODO: check POS
        sha_gamma.Update(s_square[k].elem);
        sha_gamma.Update(b_square[k].elem);
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
  sha_h.Final(h_);

  // 1.e + 1.f + 1.g
  std::vector<std::vector<FieldType>> alpha_computed;
  std::vector<FieldType> alpha_sum_computed;
  alpha_computed.resize(m);
  alpha_sum_computed.resize(m);

  for (auto mm = 0; mm < m; ++mm) {
    alpha_computed[mm].resize(N);
    alpha_sum_computed[mm] = FieldType(0);
  }

  osuCrypto::SHA1 sha_pi(sizeof(block));

  for (auto i = 0; i < N; ++i) {
    for (auto k = 0; k < m; ++k) {
      if (i != i_bar_) {
        if (i != N - 1) {
          alpha_computed[k][i] = s_[k][i] - ep_[k] * b_[k][i];

          alpha_sum_computed[k] += alpha_computed[k][i]; // 1.i

          sha_pi.Update(alpha_computed[k][i].elem);
        }
        else {
          alpha_computed[k][i] = s[k] - ep_[k] * b_[k][i];

          alpha_sum_computed[k] += alpha_computed[k][i]; // 1.i

          sha_pi.Update(alpha_computed[k][i].elem);
        }
      }
      else {
        alpha_sum_computed[k] += alpha_i_bar[k];

        sha_pi.Update(alpha_i_bar[k].elem);
      }
    }
  }
  sha_pi.Update(g_);
  sha_pi.Final(pi_);

  // 1.h + 1.i
  o_.resize(N);

  osuCrypto::SHA1 sha_psi(sizeof(block));

  for (auto i = 0; i < N; ++i) {
    o_[i] = FieldType(0);

    if (i == i_bar_) {
      o_[i] = o_i_bar;

      sha_psi.Update(o_[i].elem); // For step 2.b

      continue;
    }

    for (auto l = 0; l < n; ++l) {
      FieldType tmp(0);

      for (auto k = 0; k < m; ++k) {
        if (i != N - 1)
          tmp += a_[l][k] * s_[k][i];
        else
          tmp += a_[l][k] * s[k];
      }

      o_[i] += be_[l] * ((t_[l] / FieldType(N)) - tmp);
    }

    for (auto k = 0; k < m; ++k) {
      if (i != N - 1)
        o_[i] += ga_[k] * (s_square_[k][i] - s_[k][i]);
      else
        o_[i] += ga_[k] * (s_square[k] - s[k]);
    }

    sha_psi.Update(o_[i].elem); // For step 2.b
  }

  sha_psi.Update(w_);
  sha_psi.Final(psi_);

  // 1.j
  sigma_o = FieldType(0);

  for (auto i = 0; i < N; ++i) {
    sigma_o += o_[i];
  }

  if (sigma_o != FieldType(0)) {
    reject_ = true;

    return false;
  }

  // 1.k
  v_.resize(N);

  de_.resize(m);
  prng_seed_global_.SetSeed(seed_global_.b);

  for (auto k = 0; k < m; ++k) {
    de_[k] = FieldType(prng_seed_global_.get<block>().halves[0]);
  }

  osuCrypto::SHA1 sha_theta(sizeof(block)); // For step 1.l

  for (auto i = 0; i < N; ++i) {
    v_[i] = FieldType(0);

    if (i == i_bar_) {
      v_[i] = v_i_bar;

      sha_theta.Update(v_[i].elem); // For step 2.b

      continue;
    }

    for (auto k = 0; k < m; ++k) {
      if (i != N - 1)
        v_[i] += de_[k] * (s_square_[k][i] - alpha_sum_computed[k] * (s_[k][i] + ep_[k] * b_[k][i]) - ((ep_[k] * ep_[k]) * b_square_[k][i]));
      else
        v_[i] += de_[k] * (s_square[k] - alpha_sum_computed[k] * (s[k] + ep_[k] * b_[k][i]) - ((ep_[k] * ep_[k]) * b_square[k]));
    }

    sha_theta.Update(v_[i].elem); // For step 1.l
  }

  sha_theta.Update(u_);
  sha_theta.Final(theta_);

  // 1.m
  sigma_v = FieldType(0);

  for (auto i = 0; i < N; ++i) {
    sigma_v += v_[i];
  }

  if (sigma_v != FieldType(0)) {
    reject_ = true;

    return false;
  }

  return true;
}


}


#endif // __LZKP_SAC_VERIFIER_H_FILE___
