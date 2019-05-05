//
// Created by roee on 9/3/18.
//

#ifndef __LZKP_CAC_PROVER_LOGIC_H_FILE__
#define __LZKP_CAC_PROVER_LOGIC_H_FILE__


#include "parameters.h"
#include "cac_prover.h"

#include <thread>


namespace lzkp {


template <class FieldType>
class CacProverLogic {
public:
  CacProverLogic(const Parameters &s, FieldType **&a, FieldType *&t, FieldType *&secret, int nthreads = 1);
  ~CacProverLogic();

  void r1(block &h_gamma);
  void r3(const std::vector<uint8_t> &E, std::vector<block> &seed, std::vector<block> &omegaN, block &h_pi);
  void r5(const block &seed_ell, block &h_psi);
  void r7(const std::vector<int> &i_bar, block &seed_e_bar, std::vector<std::vector<block>> &seed_tree,
          std::vector<block> &gamma_i_bar, std::vector<std::vector<FieldType>> &alpha_i_bar, std::vector<FieldType> &o_i_bar,
          std::vector<std::vector<FieldType>> &b_square, std::vector<std::vector<FieldType>> &s);

private:
public:
  // Public known values
  FieldType **&a_;
  FieldType *&t_;

  // CacProver's secret
  FieldType *&secret_;

  const Parameters &par_;
  const int M;
  const int tau;
  const int N;
  const int n;
  const int m;

  size_t nthreads_;

  std::vector<CacProver<FieldType> *> provers_;

  std::vector<block> master_seed_;
  block h_gamma_;

  std::vector<uint8_t> E_;
  block seed_e_bar_;
  osuCrypto::PRNG prng_e_bar_;
  std::vector<block> g_;
  block h_pi_;

  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;
  block h_psi_;

  int time_eq_1;
  int tot_matrix_multiplication_time;
};

template <class FieldType>
CacProverLogic<FieldType>::CacProverLogic(const Parameters &s, FieldType **&a, FieldType *&t, FieldType *&secret, int nthreads)
    : a_(a), t_(t), secret_(secret), par_(s), M(s.M), tau(s.tau), N(s.N), n(s.n), m(s.m), nthreads_(nthreads) {
  provers_.resize(M);
}

template <class FieldType>
CacProverLogic<FieldType>::~CacProverLogic() {
  for (auto e = 0; e < M; ++e) {
    if (provers_[e])
      delete provers_[e];
  }
}

template <class FieldType>
void CacProverLogic<FieldType>::r1(block &h_gamma) {
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
            provers_[e] = new CacProver<FieldType>(par_, a_, t_, secret_);

            provers_[e]->r1(master_seed_[e]);
          }
        }, t * M / nthreads_, (t + 1) == nthreads_ ? M : (t + 1) * M / nthreads_, t));
  }
  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });

  osuCrypto::Blake2 blake_h_gamma(sizeof(block));
  std::vector<block> h_to_hash(M);
  for (auto e = 0; e < M; ++e) {
    h_to_hash[e] = provers_[e]->h_;
  }
  blake_h_gamma.Update(h_to_hash.data(), M);
  blake_h_gamma.Final(h_gamma_.bytes);

  h_gamma = h_gamma_; // Set out variable
}

template <class FieldType>
void CacProverLogic<FieldType>::r3(const std::vector<uint8_t> &E, std::vector<block> &seed, std::vector<block> &omegaN, block &h_pi) {
  E_ = E;

  // 1
  seed_e_bar_.b = osuCrypto::sysRandomSeed();
  prng_e_bar_.SetSeed(seed_e_bar_.b);

  // 2
  for (auto e = 0; e < M; ++e) {
    if (!E[e]) {
      // 2.a
      provers_[e]->g_ = prng_e_bar_.get<block>();
    }
  }

  // ** can be parallelized **
  std::vector<std::thread> threads(nthreads_);

  for(auto t = 0u; t < nthreads_; t++) {
    threads[t] = std::thread(std::bind(
        [&](const int bi, const int ei, const int t) {
          for (auto e = bi; e < ei; ++e) {
            if (!E_[e]) {
              provers_[e]->r3();
            }
          }
        }, t * M / nthreads_, (t + 1) == nthreads_ ? M : (t + 1) * M / nthreads_, t));
  }
  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });

  // 2.h
  seed.resize(tau);
  omegaN.resize(M - tau);

  osuCrypto::Blake2 blake_h_pi(sizeof(block)); // For step 3
  std::vector<block> pi_to_hash(M - tau);
  for (auto e = 0, o_id = 0, s_id = 0; e < M; ++e) {
    if (E[e]) {
      seed[s_id++] = master_seed_[e]; // For out variable
      continue;
    }

    pi_to_hash[o_id] = provers_[e]->pi_;
    omegaN[o_id++] = provers_[e]->omegaN_; // Set out variable
  }
  blake_h_pi.Update(pi_to_hash.data(), M - tau);
  blake_h_pi.Final(h_pi_);

  h_pi = h_pi_; // Set out variable
}

template <class FieldType>
void CacProverLogic<FieldType>::r5(const block &seed_ell, block &h_psi) {
  seed_ell_ = seed_ell;

  prng_seed_ell_.SetSeed(seed_ell_.b);

  for (auto e = 0; e < M; ++e) {
    if (E_[e]) {
      continue;
    }

    osuCrypto::Blake2 blake_psi(sizeof(block));

    // 1
    provers_[e]->coefficients_.resize(n + m);

    for (auto i = 0; i < n; ++i) {
      provers_[e]->coefficients_[i] = FieldType(prng_seed_ell_.get<block>().halves[0]);
    }
    for (auto i = 0; i < m; ++i) {
      provers_[e]->coefficients_[n + i] = FieldType(prng_seed_ell_.get<block>().halves[0]);
    }

    provers_[e]->w_ = prng_e_bar_.get<block>();
  }

  // ** can be parallelized **
  std::vector<std::thread> threads(nthreads_);

  for(auto t = 0u; t < nthreads_; t++) {
    threads[t] = std::thread(std::bind(
        [&](const int bi, const int ei, const int t) {
          for (auto e = bi; e < ei; ++e) {
            if (!E_[e]) {
              provers_[e]->r5();
            }
          }
        }, t * M / nthreads_, (t + 1) == nthreads_ ? M : (t + 1) * M / nthreads_, t));
  }
  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });

  osuCrypto::Blake2 blake_h_psi(sizeof(block));

  time_eq_1 = 0;
  tot_matrix_multiplication_time = 0;
  std::vector<block> psi_to_hash(M - tau);
  for (auto e = 0, e_id = 0; e < M; ++e) {
    if (!E_[e]) {
      psi_to_hash[e_id++] = provers_[e]->psi_;

      time_eq_1 += provers_[e]->time_eq_1_;
      tot_matrix_multiplication_time += provers_[e]->tot_matrix_multiplication_time_;
    }
  }
  blake_h_psi.Update(psi_to_hash.data(), M - tau);
  blake_h_psi.Final(h_psi_);

  h_psi = h_psi_;
}

template <class FieldType>
void CacProverLogic<FieldType>::r7(const std::vector<int> &i_bar, block &seed_e_bar, std::vector<std::vector<block>> &seed_tree,
                       std::vector<block> &gamma_i_bar, std::vector<std::vector<FieldType>> &alpha_i_bar, std::vector<FieldType> &o_i_bar,
                       std::vector<std::vector<FieldType>> &b_square, std::vector<std::vector<FieldType>> &s) {
  seed_tree.resize(M - tau);
  gamma_i_bar.resize(M - tau);
  alpha_i_bar.resize(M - tau);
  o_i_bar.resize(M - tau);
  b_square.resize(M - tau);
  s.resize(M - tau);

  std::vector<int> map;

  map.resize(M);
  for (auto e = 0, e_it = 0; e < M; ++e) {
    if (!E_[e]) {
      map[e] = e_it;
      provers_[e]->i_bar_ = i_bar[e_it];

      e_it++;
    }
  }

  seed_e_bar = seed_e_bar_;

  // ** can be parallelized **
  std::vector<std::thread> threads(nthreads_);

  for(auto t = 0u; t < nthreads_; t++) {
    threads[t] = std::thread(std::bind(
        [&](const int bi, const int ei, const int t) {
          for (auto e = bi; e < ei; ++e) {
            if (!E_[e]) {
              provers_[e]->r7(seed_tree[map[e]], gamma_i_bar[map[e]], alpha_i_bar[map[e]],
                              o_i_bar[map[e]], b_square[map[e]], s[map[e]]);
            }
          }
        }, t * M / nthreads_, (t + 1) == nthreads_ ? M : (t + 1) * M / nthreads_, t));
  }
  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });
}


}


#endif // __LZKP_CAC_PROVER_LOGIC_H_FILE__
