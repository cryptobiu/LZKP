//
// Created by roee on 9/3/18.
//

#ifndef __LZKP_SAC_PROVER_LOGIC_H_FILE__
#define __LZKP_SAC_PROVER_LOGIC_H_FILE__


#include "parameters.h"
#include "sac_prover.h"

#include <thread>


namespace lzkp {


template <class FieldType>
class SacProverLogic {
public:
  SacProverLogic(const Parameters &s, FieldType **&a, FieldType *&t, FieldType *&secret, bool multi_threaded = false);
  ~SacProverLogic();

  void r1(block &h_gamma);
  void r3(const block &seed_ell, block &h_pi, block &h_psi, block &h_theta);
  void r5(const std::vector<int> &i_bar, block &seed_global, std::vector<std::vector<block>> &seed_tree,
          std::vector<block> &gamma_i_bar, std::vector<std::vector<FieldType>> &alpha_i_bar, std::vector<FieldType> &o_i_bar,
          std::vector<FieldType> &v_i_bar, std::vector<std::vector<FieldType>> &b_square,
          std::vector<std::vector<FieldType>> &s, std::vector<std::vector<FieldType>> &s_square);

private:
public:
  // Public known values
  FieldType **&a_;
  FieldType *&t_;

  // CacProver's secret
  FieldType *&secret_;

  const Parameters &par_;
  const int M;
  const int N;
  const int n;
  const int m;

  size_t nthreads_;

  std::vector<SacProver<FieldType> *> provers_;

  std::vector<block> master_seed_;
  block h_gamma_;

  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;
  block seed_global_;
  osuCrypto::PRNG prng_seed_global_;
  block h_pi_, h_psi_, h_theta_;
//  std::vector<uint8_t> E_;
//  block seed_e_bar_;
//  osuCrypto::PRNG prng_e_bar_;
//  std::vector<block> g_;
//  block h_pi_;
//
//
//  block h_psi_;
};

template <class FieldType>
SacProverLogic<FieldType>::SacProverLogic(const Parameters &s, FieldType **&a, FieldType *&t, FieldType *&secret, bool multi_threaded)
    : a_(a), t_(t), secret_(secret), par_(s), M(s.M), N(s.N), n(s.n), m(s.m) {
  if (multi_threaded)
    nthreads_ = std::thread::hardware_concurrency();
  else
    nthreads_ = 1;

  provers_.resize(M);
}

template <class FieldType>
SacProverLogic<FieldType>::~SacProverLogic() {
  for (auto e = 0; e < M; ++e) {
    if (provers_[e])
      delete provers_[e];
  }
}

template <class FieldType>
void SacProverLogic<FieldType>::r1(block &h_gamma) {
  master_seed_.resize(M);

  osuCrypto::PRNG prng;

  // 1
  for (auto e = 0; e < M; ++e) {
    // 1.a
    master_seed_[e].b = osuCrypto::sysRandomSeed();
  }

  // 2 - ** can be parallelized **
  std::vector<std::thread> threads(nthreads_);

  for(auto t = 0u; t < nthreads_; t++) {
    threads[t] = std::thread(std::bind(
        [&](const int bi, const int ei, const int t) {
          for (auto e = bi; e < ei; ++e) {
            provers_[e] = new SacProver<FieldType>(par_, a_, t_, secret_);

            provers_[e]->r1(master_seed_[e]);
          }
        }, t * M / nthreads_, (t + 1) == nthreads_ ? M : (t + 1) * M / nthreads_, t));
  }
  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });

  osuCrypto::SHA1 sha_h_gamma(sizeof(block));
  for (auto e = 0; e < M; ++e) {
    sha_h_gamma.Update(provers_[e]->h_);
  }
  sha_h_gamma.Final(h_gamma_.bytes);

  h_gamma = h_gamma_; // Set out variable
}

template <class FieldType>
void SacProverLogic<FieldType>::r3(const block &seed_ell, block &h_pi, block &h_psi, block &h_theta) {
  seed_ell_ = seed_ell;
  prng_seed_ell_.SetSeed(seed_ell_.b);

  // 2
  seed_global_.b = osuCrypto::sysRandomSeed();
  prng_seed_global_.SetSeed(seed_global_.b);

  for (auto e = 0; e < M; ++e) {
    provers_[e]->seed_ell_ = prng_seed_ell_.get<block>();
    provers_[e]->g_ = prng_seed_global_.get<block>();
    provers_[e]->w_ = prng_seed_global_.get<block>();
    provers_[e]->u_ = prng_seed_global_.get<block>();

    provers_[e]->seed_global_ = prng_seed_global_.get<block>();
  }

  // ** can be parallelized **
  std::vector<std::thread> threads(nthreads_);

  for(auto t = 0u; t < nthreads_; t++) {
    threads[t] = std::thread(std::bind(
        [&](const int bi, const int ei, const int t) {
          for (auto e = bi; e < ei; ++e) {
            provers_[e]->r3();
          }
        }, t * M / nthreads_, (t + 1) == nthreads_ ? M : (t + 1) * M / nthreads_, t));
  }
  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });

  // 4 + 5 + 6
  osuCrypto::SHA1 sha_h_pi(sizeof(block)); // For step 4
  osuCrypto::SHA1 sha_h_psi(sizeof(block)); // For step 5
  osuCrypto::SHA1 sha_h_theta(sizeof(block)); // For step 6

  for (auto e = 0; e < M; ++e) {
    sha_h_pi.Update(provers_[e]->pi_);
    sha_h_psi.Update(provers_[e]->psi_);
    sha_h_theta.Update(provers_[e]->theta_);
  }

  sha_h_pi.Final(h_pi_); // 4
  sha_h_psi.Final(h_psi_); // 5
  sha_h_theta.Final(h_theta_); // 6

  h_pi = h_pi_;
  h_psi = h_psi_;
  h_theta = h_theta_;
}

template <class FieldType>
void SacProverLogic<FieldType>::r5(const std::vector<int> &i_bar, block &seed_global, std::vector<std::vector<block>> &seed_tree,
               std::vector<block> &gamma_i_bar, std::vector<std::vector<FieldType>> &alpha_i_bar, std::vector<FieldType> &o_i_bar,
               std::vector<FieldType> &v_i_bar, std::vector<std::vector<FieldType>> &b_square,
               std::vector<std::vector<FieldType>> &s, std::vector<std::vector<FieldType>> &s_square) {
  seed_tree.resize(M);
  gamma_i_bar.resize(M);
  alpha_i_bar.resize(M);
  o_i_bar.resize(M);
  v_i_bar.resize(M);
  b_square.resize(M);
  s.resize(M);
  s_square.resize(M);

  seed_global = seed_global_;

  for (auto e = 0; e < M; ++e) {
    provers_[e]->i_bar_ = i_bar[e];
  }
  // ** can be parallelized **
  std::vector<std::thread> threads(nthreads_);

  for(auto t = 0u; t < nthreads_; t++) {
    threads[t] = std::thread(std::bind(
        [&](const int bi, const int ei, const int t) {
          for (auto e = bi; e < ei; ++e) {
              provers_[e]->r5(seed_tree[e], gamma_i_bar[e], alpha_i_bar[e],
                              o_i_bar[e], v_i_bar[e], b_square[e], s[e], s_square[e]);
          }
        }, t * M / nthreads_, (t + 1) == nthreads_ ? M : (t + 1) * M / nthreads_, t));
  }
  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });
}


}


#endif // __LZKP_SAC_PROVER_LOGIC_H_FILE__