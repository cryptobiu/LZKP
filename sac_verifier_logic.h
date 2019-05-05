//
// Created by roee on 9/3/18.
//

#ifndef __LZKP_SAC_VERIFIER_LOGIC_H_FILE__
#define __LZKP_SAC_VERIFIER_LOGIC_H_FILE__


#include "parameters.h"
#include "sac_verifier.h"

#include <thread>


namespace lzkp {


template <class FieldType>
class SacVerifierLogic {
public:
  SacVerifierLogic(const Parameters &s, FieldType **&a, FieldType *&t, int nthreads = 1);
  ~SacVerifierLogic();

  void r2(const block &h_gamma, block &seed_ell);
  void r4(const block &h_pi, const block &h_psi, const block &h_theta, std::vector<int> &i_bar);
  bool r6(const block &seed_global, const std::vector<std::vector<block>> &seed_tree, const std::vector<block> &gamma_i_bar,
          const std::vector<std::vector<FieldType>> &alpha_i_bar, const std::vector<FieldType> &o_i_bar, const std::vector<FieldType> &v_i_bar,
          const std::vector<std::vector<FieldType>> &b_square, const std::vector<std::vector<FieldType>> &s, const std::vector<std::vector<FieldType>> &s_square);

private:
public:
  // Public known values
  FieldType **&a_;
  FieldType *&t_;

  const Parameters &par_;
  const int M;
  const int N;
  const int n;
  const int m;

  size_t nthreads_;

  std::vector<SacVerifier<FieldType> *> verifiers_;

  block h_gamma_;
  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;

  block h_pi_;
  block h_psi_;
  block h_theta_;

  block seed_global_;
  osuCrypto::PRNG prng_seed_global_;

  int time_eq_1;
  long tot_matrix_multiplication_time;
};


template <class FieldType>
SacVerifierLogic<FieldType>::SacVerifierLogic(const Parameters &s, FieldType **&a, FieldType *&t, int nthreads)
    : a_(a), t_(t), par_(s), M(s.M), N(s.N), n(s.n), m(s.m), nthreads_(nthreads) {
  verifiers_.resize(M);
}

template <class FieldType>
SacVerifierLogic<FieldType>::~SacVerifierLogic() {
  for (auto e = 0; e < M; ++e) {
    if (verifiers_[e])
      delete verifiers_[e];
  }
}

template <class FieldType>
void SacVerifierLogic<FieldType>::r2(const block &h_gamma, block &seed_ell) {
  h_gamma_ = h_gamma;

  seed_ell_.b = osuCrypto::sysRandomSeed();
  prng_seed_ell_.SetSeed(seed_ell_.b);

  for (auto e = 0; e < M; ++e) {
    verifiers_[e] = new SacVerifier<FieldType>(par_, a_, t_);

    verifiers_[e]->seed_ell_ = prng_seed_ell_.get<block>();
  }

  // ** can be parallelized **
  std::vector<std::thread> threads(nthreads_);

  for(auto t = 0u; t < nthreads_; t++) {
    threads[t] = std::thread(std::bind(
        [&](const int bi, const int ei, const int t) {
          for (auto e = bi; e < ei; ++e) {
            verifiers_[e]->r2();
          }
        }, t * M / nthreads_, (t + 1) == nthreads_ ? M : (t + 1) * M / nthreads_, t));
  }
  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });

  seed_ell = seed_ell_;
}

template <class FieldType>
void SacVerifierLogic<FieldType>::r4(const block &h_pi, const block &h_psi, const block &h_theta, std::vector<int> &i_bar) {
  h_pi_ = h_pi;
  h_psi_ = h_psi;
  h_theta_ = h_theta;

  i_bar.resize(M);

  osuCrypto::PRNG prng;
  prng.SetSeed(osuCrypto::sysRandomSeed());

  for (auto e = 0; e < M; ++e) {
    verifiers_[e]->i_bar_ = prng.get<osuCrypto::u32>() % N;

    i_bar[e] = verifiers_[e]->i_bar_; // Set out variable
  }
}

template <class FieldType>
bool SacVerifierLogic<FieldType>::r6(const block &seed_global, const std::vector<std::vector<block>> &seed_tree, const std::vector<block> &gamma_i_bar,
                             const std::vector<std::vector<FieldType>> &alpha_i_bar, const std::vector<FieldType> &o_i_bar, const std::vector<FieldType> &v_i_bar,
                             const std::vector<std::vector<FieldType>> &b_square, const std::vector<std::vector<FieldType>> &s, const std::vector<std::vector<FieldType>> &s_square) {
  seed_global_ = seed_global;
  prng_seed_global_.SetSeed(seed_global_.b);

  // Generate g_[e] = for step 1.g
  for (auto e = 0; e < M; ++e) {
    verifiers_[e]->g_ = prng_seed_global_.get<block>();
    verifiers_[e]->w_ = prng_seed_global_.get<block>();
    verifiers_[e]->u_ = prng_seed_global_.get<block>();

    verifiers_[e]->seed_global_ = prng_seed_global_.get<block>();
  }

  // ** can be parallelized **
  std::vector<std::thread> threads(nthreads_);

  for(auto t = 0u; t < nthreads_; t++) {
    threads[t] = std::thread(std::bind(
        [&](const int bi, const int ei, const int t) {
          for (auto e = bi; e < ei; ++e) {
            verifiers_[e]->r6(seed_tree[e], gamma_i_bar[e], alpha_i_bar[e], o_i_bar[e], v_i_bar[e],
                              b_square[e], s[e], s_square[e]);
          }
        }, t * M / nthreads_, (t + 1) == nthreads_ ? M : (t + 1) * M / nthreads_, t));
  }
  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });

  for (auto e = 0; e < M; ++e) {
    if (verifiers_[e]->reject_)
      return false;
  }

  // 2 + 3 + 4 + 5
  osuCrypto::Blake2 blake_h_gamma(sizeof(block));
  osuCrypto::Blake2 blake_h_pi(sizeof(block));
  osuCrypto::Blake2 blake_h_psi(sizeof(block));
  osuCrypto::Blake2 blake_h_theta(sizeof(block));

  std::vector<block> h_to_hash(M);
  std::vector<block> pi_to_hash(M);
  std::vector<block> psi_to_hash(M);
  std::vector<block> theta_to_hash(M);

  time_eq_1 = 0;
  tot_matrix_multiplication_time = 0;
  for (auto e = 0; e < M; ++e) {
    h_to_hash[e] = verifiers_[e]->h_; // For step 2
    pi_to_hash[e] = verifiers_[e]->pi_; // For step 3
    psi_to_hash[e] = verifiers_[e]->psi_; // For step 4
    theta_to_hash[e] = verifiers_[e]->theta_; // For step 5

    time_eq_1 += verifiers_[e]->time_eq_1_;
    tot_matrix_multiplication_time += verifiers_[e]->tot_matrix_multiplication_time_ / 1000000;
  }

  // 2
  block gamma;
  blake_h_gamma.Update(h_to_hash.data(), M);
  blake_h_gamma.Final(gamma);
  if (!eq(h_gamma_.b, gamma.b))
    return false;

  // 3
  block pi;
  blake_h_pi.Update(pi_to_hash.data(), M);
  blake_h_pi.Final(pi);
  if (!eq(pi.b, h_pi_.b))
    return false;

  // 4
  block psi;
  blake_h_psi.Update(psi_to_hash.data(), M);
  blake_h_psi.Final(psi);
  if (!eq(psi.b, h_psi_.b))
    return false;

  // 5
  block theta;
  blake_h_theta.Update(theta_to_hash.data(), M);
  blake_h_theta.Final(theta);
  if (!eq(theta.b, h_theta_.b))
    return false;

  return true;
}


}


#endif // __LZKP_SAC_VERIFIER_LOGIC_H_FILE__
