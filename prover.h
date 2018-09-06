#ifndef __LZKP_PROVER_H_FILE___
#define __LZKP_PROVER_H_FILE___


#include <stack>
#include <cryptoTools/Crypto/sha1.h>

#include "seedtree.h"
#include "settings.h"

//#include <NTL/ZZ_p.h>


namespace lzkp {


template <class FieldType>
class Prover {
public:
  Prover(const Settings &s, const std::vector<std::vector<FieldType>> &a, const std::vector<FieldType> &t, const std::vector<FieldType> &secret);
  ~Prover();

  void r1(const block &master_seed);
  void r3(); // Local variable g_ must be set before calling this method
  void r5(); // Local variable coefficients_ must be set before calling this method
  void r7(std::vector<block> &seed_tree, block &gamma_i_bar, std::vector<FieldType> &alpha_i_bar,
          FieldType &o_i_bar, std::vector<FieldType> &b_square, std::vector<FieldType> &s);
             // Local variable i_bar_ must be set before calling this method

private:
public:
  // Public known values
  const std::vector<std::vector<FieldType>> &a_;
  const std::vector<FieldType> &t_;

  // Prover's secret
  const std::vector<FieldType> secret_;

  const int N;
//  const uint64_t q;
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
};

template <class FieldType>
Prover<FieldType>::Prover(const Settings &s, const std::vector<std::vector<FieldType>> &a, const std::vector<FieldType> &t, const std::vector<FieldType> &secret)
    : a_(a), t_(t), secret_(secret), N(s.N), n(s.n), m(s.m) {
}

template <class FieldType>
Prover<FieldType>::~Prover() {
}

template <class FieldType>
void Prover<FieldType>::r1(const block &master_seed) {
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
      b_square_[k][N - 1] += b_[k][i]; // Accumulate b_e,k
    }
    b_square_[k][N - 1] *= b_square_[k][N - 1]; // Calculate b_e,k^2
    for (auto i = 0; i < N - 1; ++i) {
      b_square_[k][N - 1] -= b_square_[k][i]; // Subtract b^2_e,k,i
    }
  }

  // 1.e
  gamma_.resize(N);
  osuCrypto::SHA1 sha_gamma(sizeof(block));

  for (auto i = 0; i < N; ++i) {
    block blk = seed_tree_.getSeed(i);
    sha_gamma.Reset(); // Calculate com(state_e,i , r_e,i) == com(seed_e,i , r_e,i)
    sha_gamma.Update(blk); // Hash seed_e,i
    if (i == N - 1) {
//      unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT
      for (auto k = 0; k < m; ++k) { // TODO: optimize
//        NTL::BytesFromZZ(buf, b_square_[k][N - 1].elem._ZZ_p__rep, NTL::NumBytes(b_square_[k][N - 1].elem._ZZ_p__rep));
//        sha_gamma.Update(buf, NTL::NumBytes(b_square_[k][N - 1].elem._ZZ_p__rep));
          sha_gamma.Update(b_square_[k][N - 1].elem); // TODO: check POS
      }
    }
    sha_gamma.Update(r_[i]); // Hash r_e,i
    sha_gamma.Final(gamma_[i]);
//    std::cout << "Gamma " << i << " " << gamma_[i].halves[0] << " " << gamma_[i].halves[1] << std::endl;
  }

  // 1.f
  osuCrypto::SHA1 sha_h(sizeof(block));
  for (auto i = 0; i < N; ++i) {
    sha_h.Update(gamma_[i]);
  }
  sha_h.Final(h_); // mb need to zero it first
}

template <class FieldType>
void Prover<FieldType>::r3() {
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
      //std::cout << "alpha " << k << " " << i << " " << alpha_[k][i] << std::endl;
    }
    //std::cout << "alpha_sum " << k << " " << alpha_sum_[k] << std::endl;
  }

  // 2.f
  osuCrypto::SHA1 sha_omegaN(sizeof(block));
//  unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT
  for (auto i = 0; i < m; ++i) {
//    NTL::BytesFromZZ(buf, s_[i][N - 1].elem._ZZ_p__rep, NTL::NumBytes(s_[i][N - 1].elem._ZZ_p__rep));
//    sha_omegaN.Update(buf, NTL::NumBytes(s_[i][N - 1].elem._ZZ_p__rep));
      sha_omegaN.Update(s_[i][N - 1].elem);
  }
  sha_omegaN.Update(gN_);
  sha_omegaN.Final(omegaN_);

  // 2.g
  osuCrypto::SHA1 sha_pi(sizeof(block));

  for (auto mm = 0; mm < m; ++mm) {
    for (auto nn = 0; nn < N; ++nn) {
//      NTL::BytesFromZZ(buf, alpha_[mm][nn].elem._ZZ_p__rep, NTL::NumBytes(alpha_[mm][nn].elem._ZZ_p__rep));
//      sha_pi.Update(buf, NTL::NumBytes(alpha_[mm][nn].elem._ZZ_p__rep));
      sha_pi.Update(alpha_[mm][nn].elem);
    }
  }
  sha_pi.Update(g_);
  sha_pi.Final(pi_);

  //std::cout << "PI E " << e << " " << pi.halves[0] << pi.halves[1] << std::endl;
}

template <class FieldType>
void Prover<FieldType>::r5() {
  o_.resize(N);

  // 2.a
  osuCrypto::SHA1 sha_psi(sizeof(block));
//  unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT

  for (auto i = 0; i < N; ++i) {
    o_[i] = FieldType(0);

    for (auto l = 0; l < n; ++l) {
      FieldType tmp(0);

      for (auto k = 0; k < m; ++k) {
        tmp += a_[l][k] * s_[k][i];
      }

      o_[i] += coefficients_[l] * ((t_[l] / FieldType(N)) - tmp);
//      o_[i] += coefficients_[l] * (t_[l] - tmp);
    }

    for (auto k = 0; k < m; ++k) {
      //o_[e_it][i] += coefficients_[e_it][n + k] *
//                       (alpha_sum_[e][k] * (s_[e][k][i] + b_[e][k][i]) + b_square_[e][k][i] - s_[e][k][i]);
      o_[i] += coefficients_[n + k] *
               (alpha_sum_[k] * (s_[k][i] + b_[k][i]) + b_square_[k][i] - s_[k][i]);
    }

    // 2.b
//    NTL::BytesFromZZ(buf, o_[i].elem._ZZ_p__rep, NTL::NumBytes(o_[i].elem._ZZ_p__rep));
//    sha_psi.Update(buf, NTL::NumBytes(o_[i].elem._ZZ_p__rep)); // For step 2.b
    sha_psi.Update(o_[i].elem); // For step 2.b

//    std::cout << "O " << e << " " << i << " " << o_[e_it][i] << std::endl;
  }

  sha_psi.Update(w_);
  sha_psi.Final(psi_);
}

template <class FieldType>
void Prover<FieldType>::r7(std::vector<block> &seed_tree, block &gamma_i_bar, std::vector<FieldType> &alpha_i_bar,
                FieldType &o_i_bar, std::vector<FieldType> &b_square, std::vector<FieldType> &s) {
  seed_tree.clear();

  std::stack<std::pair<int, std::pair<int, int>>> seed_stack;

  seed_stack.push(std::make_pair(0, std::make_pair(0, N)));

//  std::cout << i_bar_ << "\t";
  while (!seed_stack.empty()) {
    auto item = seed_stack.top();
    seed_stack.pop();

    // Check for leaf
    if (item.second.second - item.second.first == 1) {
      if (item.second.first != i_bar_) {
//        std::cout << "pushing leaf " << item.second.first << " " << seed_tree_.getSeed(item.second.first).halves[0] << std::endl;
        seed_tree.push_back(seed_tree_.getSeed(item.second.first));
      }
      continue;
    }

    if (i_bar_ >= item.second.first && i_bar_ < item.second.second) { // seed_e_i_bar[e] is a descendant of the current node
      seed_stack.push(std::make_pair(item.first * 2 + 2, std::make_pair((item.second.first + item.second.second) / 2, item.second.second))); // Add right child
      seed_stack.push(std::make_pair(item.first * 2 + 1, std::make_pair(item.second.first, (item.second.first + item.second.second) / 2))); // Add left child
    }
    else {
//      std::cout << "pushing seed " << item.first << " " << seed_tree_[item.first].halves[0] << std::endl;
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
#endif