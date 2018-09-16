#ifndef __LZKP_CAC_VERIFIER_H_FILE___
#define __LZKP_CAC_VERIFIER_H_FILE___


#include <stack>
#include <queue>
#include <cryptoTools/Crypto/Blake2.h>
#include <cryptoTools/Crypto/PRNG.h>

#include "seedtree.h"
#include "parameters.h"


namespace lzkp {


template <class FieldType>
class CacVerifier {
public:
  CacVerifier(const Parameters &s, FieldType **&a, FieldType *&t);
  ~CacVerifier();

  void r4(); // Local variable seed_ must be set before calling this method
  bool r8(const std::vector<block> &seed_tree, const block &gamma_i_bar,
          const std::vector<FieldType> &alpha_i_bar, const FieldType &o_i_bar, const std::vector<FieldType> &b_square,
          const std::vector<FieldType> &s); // Local variables omegaN_ (from round 4), g_, w_,  must be set before calling this method

private:
public:
  // Public known values
  FieldType **&a_;
  FieldType *&t_;

  const int N;
  const int n;
  const int m;

  // Local values
  block seed_;
  block omegaN_;
  SeedTree seed_tree_;
  std::vector<block> r_;
  std::vector<std::vector<FieldType>> b_;
  std::vector<std::vector<FieldType>> b_square_;
  std::vector<block> gamma_;
  block h_;
  std::vector<FieldType> coefficients_;

  int i_bar_;

  bool reject_;
  std::vector<block> partial_seeds_;
  block pi_;
  block g_;
  std::vector<FieldType> o_;
  FieldType sigma_o;
  block w_;
  block psi_;
};

template <class FieldType>
CacVerifier<FieldType>::CacVerifier(const Parameters &s, FieldType **&a, FieldType *&t)
    : a_(a), t_(t), N(s.N), n(s.n), m(s.m) {
}

template <class FieldType>
CacVerifier<FieldType>::~CacVerifier() {
}

template <class FieldType>
void CacVerifier<FieldType>::r4() {
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
      b_[mm][i]        = FieldType(seed_tree_.getBlock(i).halves[0]); // NEED TO FIND A WAY TO USE ALL 128 BITS
      b_square_[mm][i] = FieldType(seed_tree_.getBlock(i).halves[0]);
    }
  }

  // *1* - 1.c
  r_[N - 1] = seed_tree_.getBlock(N - 1);

  for (auto mm = 0; mm < m; ++mm) {
    b_[mm][N - 1] = FieldType(seed_tree_.getBlock(N - 1).halves[0]);
  }

  // *1* - 1.d
  for (auto k = 0; k < m; ++k) {
    b_square_[k][N - 1] = FieldType(0);
    for (auto i = 0; i < N; ++i) {
      b_square_[k][N - 1] += b_[k][i];
    }
    b_square_[k][N - 1] *= b_square_[k][N - 1];
    for (auto i = 0; i < N - 1; ++i) {
      b_square_[k][N - 1] -= b_square_[k][i];
    }
  }

  // *1* - 1.e
  osuCrypto::Blake2 blake_gamma(sizeof(block));
  for (auto i = 0; i < N - 1; ++i) {
    block blk = seed_tree_.getSeed(i);
    blake_gamma.Reset();
    blake_gamma.Update(blk);
    blake_gamma.Update(r_[i]);
    blake_gamma.Final(gamma_[i]);
  }
  block blk = seed_tree_.getSeed(N - 1);
  blake_gamma.Reset();
  blake_gamma.Update(blk);
  for (auto i = 0; i < m; ++i) {
    blake_gamma.Update(b_square_[i][N - 1].elem);
  }
  blake_gamma.Update(r_[N - 1]);
  blake_gamma.Final(gamma_[N - 1]);

  // *1* - 1.f
  osuCrypto::Blake2 blake_h(sizeof(block));
  for (auto i = 0; i < N; ++i) {
    blake_h.Update(gamma_[i]);
  }
  blake_h.Final(h_); // mb need to zero it first
}

template <class FieldType>
bool CacVerifier<FieldType>::r8(const std::vector<block> &seed_tree, const block &gamma_i_bar,
                  const std::vector<FieldType> &alpha_i_bar, const FieldType &o_i_bar,
                  const std::vector<FieldType> &b_square, const std::vector<FieldType> &s) {
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
  gamma_.resize(N);
  r_.resize(N);
  b_.resize(m);
  b_square_.resize(m);

  for (auto mm = 0; mm < m; ++mm){
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
      b_[mm][i]        = FieldType(prng[i].get<block>().halves[0]); // NEED TO FIND A WAY TO USE ALL 128 BITS
      if (i != N - 1)
        b_square_[mm][i] = FieldType(prng[i].get<block>().halves[0]); // We ignore b_square_[e][k][N-1]
    }
  }

  // 1.c
  osuCrypto::Blake2 blake_gamma(sizeof(block));
  for (auto i = 0; i < N; ++i) {
    if (i == i_bar_)
      continue;

    block blk = partial_seeds_[i];
    blake_gamma.Reset();
    blake_gamma.Update(blk);
    if (i == N - 1) {
      for (auto k = 0; k < m; ++k) {
        blake_gamma.Update(b_square[k].elem);
      }
    }
    blake_gamma.Update(r_[i]);
    blake_gamma.Final(gamma_[i]);

  }

  // 1.d
  osuCrypto::Blake2 blake_h(sizeof(block));
  for (auto i = 0; i < N; ++i) {
    if (i != i_bar_)
      blake_h.Update(gamma_[i]);
    else
      blake_h.Update(gamma_i_bar);
  }
  blake_h.Final(h_);

  // 1.e
  if (i_bar_ != N - 1) {
    block gN_e, blk;

    gN_e = prng[N - 1].get<block>();
    osuCrypto::Blake2 blake_omega(sizeof(block));
    for (auto i = 0; i < m; ++i) {
      blake_omega.Update(s[i].elem);
    }
    blake_omega.Update(gN_e);
    blake_omega.Final(blk);

    if (!eq(blk.b, omegaN_.b)) {
      reject_ = true;

      return false;
    }
  }

  // 1.f + 1.g + 1.i
  std::vector<std::vector<FieldType>> s_computed;
  std::vector<std::vector<FieldType>> alpha_computed;
  std::vector<FieldType> alpha_sum_computed;
  s_computed.resize(m);
  alpha_computed.resize(m);
  alpha_sum_computed.resize(m);

  for (auto mm = 0; mm < m; ++mm) {
    s_computed[mm].resize(N);
    alpha_computed[mm].resize(N);
    alpha_sum_computed[mm] = FieldType(0);
  }

  // 1.f
  for (auto i = 0; i < N - 1; ++i) {
    if (i == i_bar_)
      continue;

    for (auto k = 0; k < m; ++k)
      s_computed[k][i] = FieldType(prng[i].get<block>().halves[0]); // PROBABLY NEED TO CHANGE THAT
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
  osuCrypto::Blake2 blake_pi(sizeof(block));

  for (auto mm = 0; mm < m; ++mm) {
    for (auto nn = 0; nn < N; ++nn) {
      if (nn != i_bar_) {
        blake_pi.Update(alpha_computed[mm][nn].elem);
      }
      else {
        blake_pi.Update(alpha_i_bar[mm].elem);
      }
    }
  }
  blake_pi.Update(g_);
  blake_pi.Final(pi_);

  // 1.j
  osuCrypto::Blake2 blake_psi(sizeof(block));

  o_.resize(N);
  for (auto i = 0; i < N; ++i) {
    o_[i] = FieldType(0);

    if (i == i_bar_) {
      o_[i] = o_i_bar;

      continue;
    }

    for (auto l = 0; l < n; ++l) {
      FieldType tmp(0);

      for (auto k = 0; k < m; ++k) {
        if (i != N - 1)
          tmp += a_[l][k] * s_computed[k][i];
        else
          tmp += a_[l][k] * s[k];
      }

      o_[i] += coefficients_[l] * ((t_[l] / FieldType(N)) - tmp);
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
  }

  // 1.k
  sigma_o = FieldType(0);

  for (auto i = 0; i < N; ++i) {
    blake_psi.Update(o_[i].elem);

    sigma_o += o_[i]; // for step 1.k
  }
  blake_psi.Update(w_);
  blake_psi.Final(psi_);

  // 1.l
  if (sigma_o != FieldType(0)) {
    reject_ = true;

    return false;
  }

  return true;
}


}


#endif // __LZKP_CAC_VERIFIER_H_FILE___
