#ifndef __LZKP_CAC_PROVER_H_FILE___
#define __LZKP_CAC_PROVER_H_FILE___


#include <stack>
#include <cryptoTools/Crypto/Blake2.h>
#include <chrono>

#include "seedtree.h"
#include "parameters.h"


namespace lzkp {


template <class FieldType>
class CacProver {
public:
  CacProver(const Parameters &s, FieldType **&a, FieldType *&t, FieldType *&secret);
  ~CacProver();

  void r1(const block &master_seed);
  void r3(); // Local variable g_ must be set before calling this method
  void r5(); // Local variable coefficients_ must be set before calling this method
  void r7(std::vector<block> &seed_tree, block &gamma_i_bar, std::vector<FieldType> &alpha_i_bar,
          FieldType &o_i_bar, std::vector<FieldType> &b_square, std::vector<FieldType> &s);
             // Local variable i_bar_ must be set before calling this method

private:
public:
  // Public known values
  FieldType **&a_;
  FieldType *&t_;

  // CacProver's secret
  FieldType *&secret_;

  const int N;
  const int n;
  const int m;

  // Local values
  block master_seed_;
  SeedTree seed_tree_;
  std::vector<block> r_;
  std::vector<std::vector<FieldType>> b_;
  std::vector<std::vector<FieldType>> b_square_;
  std::vector<block> gamma_;
  block h_;

  block g_;
  std::vector<std::vector<FieldType>> s_;
  std::vector<std::vector<FieldType>> alpha_;
  std::vector<FieldType> alpha_sum_;
  block gN_; // g_{e,N}
  block omegaN_;
  block pi_;

  std::vector<FieldType> coefficients_;
  std::vector<FieldType> o_;
  block w_;
  block psi_;

  int i_bar_;

  int time_eq_1_;
  long tot_matrix_multiplication_time_;
};

template <class FieldType>
CacProver<FieldType>::CacProver(const Parameters &s, FieldType **&a, FieldType *&t, FieldType *&secret)
    : a_(a), t_(t), secret_(secret), N(s.N), n(s.n), m(s.m) {
}

template <class FieldType>
CacProver<FieldType>::~CacProver() {
}

template <class FieldType>
void CacProver<FieldType>::r1(const block &master_seed) {
  master_seed_ = master_seed;

  // 1.a
  seed_tree_.resize(N);
  seed_tree_.generate(master_seed_);

  // 1.b
  r_.resize(N);
  b_.resize(m);
  b_square_.resize(m);

  for (auto mm = 0; mm < m; ++mm){
    b_[mm].resize(N);
    b_square_[mm].resize(N);
  }

  for (auto i = 0; i < N - 1; ++i) {
    r_[i] = seed_tree_.getBlock(i);

    for (auto mm = 0; mm < m; ++mm) {
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
      b_square_[k][N - 1] -= b_square_[k][i]; // Subtract b^2_k,i
    }
  }

  // 1.e
  gamma_.resize(N);
  osuCrypto::Blake2 blake_gamma(sizeof(block));

  for (auto i = 0; i < N; ++i) {
    block blk = seed_tree_.getSeed(i);
    blake_gamma.Reset(); // Calculate com(state_e,i , r_e,i) == com(seed_e,i , r_e,i)
    blake_gamma.Update(blk); // Hash seed_e,i
    if (i == N - 1) {
      std::vector<FieldType> b_square_to_hash(m);
      for (auto k = 0; k < m; ++k) { // TODO: optimize
        b_square_to_hash[k] = b_square_[k][N - 1];
      }
      blake_gamma.Update((decltype(b_square_to_hash[0].elem)*)b_square_to_hash.data(), m); // TODO: make POS
    }
    blake_gamma.Update(r_[i]); // Hash r_e,i
    blake_gamma.Final(gamma_[i]);
  }

  // 1.f
  osuCrypto::Blake2 blake_h(sizeof(block));
  blake_h.Update(gamma_.data(), N);
  blake_h.Final(h_); // mb need to zero it first
}

template <class FieldType>
void CacProver<FieldType>::r3() {
  s_.resize(m);
  alpha_.resize(m);
  alpha_sum_.resize(m);
  for (auto mm = 0; mm < m; ++mm) {
    s_[mm].resize(N);
    alpha_[mm].resize(N);
  }

  // 2.b
  for (auto i = 0; i < N - 1; ++i) {
    for (auto mm = 0; mm < m; ++mm) {
      s_[mm][i] = FieldType(seed_tree_.getBlock(i).halves[0]); // PROBABLY NEED TO CHANGE THAT
    }
  }

  // 2.c
  gN_ = seed_tree_.getBlock(N - 1);

  // 2.d
  for (auto k = 0; k < m; ++k) {
    s_[k][N - 1] = secret_[k];

    for (auto i = 0; i < N - 1; ++i) {
      s_[k][N - 1] -= s_[k][i];
    }
  }

  // 2.e + 2.h
  for (auto k = 0; k < m; ++k) {
    alpha_sum_[k] = FieldType(0);

    for (auto i = 0; i < N; ++i) {
      alpha_[k][i] = s_[k][i] - b_[k][i];
      alpha_sum_[k] += alpha_[k][i];
    }
  }

  // 2.f
  osuCrypto::Blake2 blake_omegaN(sizeof(block));
  std::vector<FieldType> s_to_hash(m);
  for (auto i = 0; i < m; ++i) {
    s_to_hash[i] = s_[i][N - 1];
  }
  blake_omegaN.Update((decltype(s_to_hash[0].elem)*)s_to_hash.data(), m); // TODO: make POS
  blake_omegaN.Update(gN_);
  blake_omegaN.Final(omegaN_);

  // 2.g
  osuCrypto::Blake2 blake_pi(sizeof(block));

  for (auto mm = 0; mm < m; ++mm) {
    blake_pi.Update((decltype(alpha_[mm][0].elem)*)alpha_[mm].data(), N);
  }
  blake_pi.Update(g_);
  blake_pi.Final(pi_);
}

template <class FieldType>
void CacProver<FieldType>::r5() {
  o_.resize(N);

  // 2.a

  auto eq_1_clock = std::chrono::high_resolution_clock::now();
  tot_matrix_multiplication_time_ = 0;

  for (auto i = 0; i < N; ++i) {
    o_[i] = FieldType(0);

    for (auto l = 0; l < n; ++l) {
      auto matrix_multiplication_start_time = std::chrono::high_resolution_clock::now();

      FieldType prod = FieldType::dotProdct(a_, s_, l, i, m);

      auto matrix_multiplication_end_time = std::chrono::high_resolution_clock::now();
      auto dur = std::chrono::duration_cast<std::chrono::nanoseconds>(matrix_multiplication_end_time - matrix_multiplication_start_time);
      tot_matrix_multiplication_time_ += dur.count();

      o_[i] += coefficients_[l] * ((t_[l] / FieldType(N)) - prod);
    }

    for (auto k = 0; k < m; ++k) {
      o_[i] += coefficients_[n + k] *
               (alpha_sum_[k] * (s_[k][i] + b_[k][i]) + b_square_[k][i] - s_[k][i]);
    }
  }

  time_eq_1_ = std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::high_resolution_clock::now() - eq_1_clock).count();

  // 2.b
  osuCrypto::Blake2 blake_psi(sizeof(block));

  blake_psi.Update((decltype(o_[0].elem)*)o_.data(), N);
  blake_psi.Update(w_);
  blake_psi.Final(psi_);
}

template <class FieldType>
void CacProver<FieldType>::r7(std::vector<block> &seed_tree, block &gamma_i_bar, std::vector<FieldType> &alpha_i_bar,
                FieldType &o_i_bar, std::vector<FieldType> &b_square, std::vector<FieldType> &s) {
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

  if (i_bar_ != N - 1) {
    b_square.resize(m);
    s.resize(m);

    for (auto k = 0; k < m; ++k) {
      b_square[k] = b_square_[k][N - 1];
      s[k] = s_[k][N - 1];
    }
  }
}


}
#endif // __LZKP_CAC_PROVER_H_FILE___