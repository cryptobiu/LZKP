//
// Created by roee on 9/3/18.
//

#ifndef LZKP_PROVER_WRAPPER_H
#define LZKP_PROVER_WRAPPER_H


#include "settings.h"
#include "prover.h"
#include "Mersenne.h"


namespace lzkp {


template <class FieldType>
class ProverWrapper {
public:
  ProverWrapper(const Settings &s, const std::vector<std::vector<FieldType>> &a, const std::vector<FieldType> &t, const std::vector<FieldType> &secret);
  ~ProverWrapper();

  void r1(block &h_gamma);
  void r3(const std::vector<bool> &E, std::vector<block> &seed, std::vector<block> &omegaN, block &h_pi);
  void r5(const block &seed_ell, block &h_psi);
  void r7(const std::vector<int> &i_bar, block &seed_e_bar, std::vector<std::vector<block>> &seed_tree,
          std::vector<block> &gamma_i_bar, std::vector<std::vector<FieldType>> &alpha_i_bar, std::vector<FieldType> &o_i_bar,
          std::vector<std::vector<FieldType>> &b_square, std::vector<std::vector<FieldType>> &s);

private:
public:
  // Public known values
  const std::vector<std::vector<FieldType>> &a_;
  const std::vector<FieldType> &t_;

  // Prover's secret
  const std::vector<FieldType> secret_;

  const Settings &set_;
  const int M;
  const int tau;
  const int N;
//  const uint64_t q;
  const int n;
  const int m;

  std::vector<Prover<FieldType> *> provers_;

  std::vector<block> master_seed_;
  block h_gamma_;

  std::vector<bool> E_;
  block seed_e_bar_;
  osuCrypto::PRNG prng_e_bar_;
  std::vector<block> g_;
  block h_pi_;

  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;
  block h_psi_;
};

template <class FieldType>
ProverWrapper<FieldType>::ProverWrapper(const Settings &s, const std::vector<std::vector<FieldType>> &a, const std::vector<FieldType> &t, const std::vector<FieldType> &secret)
    : a_(a), t_(t), secret_(secret), set_(s), M(s.M), tau(s.tau), N(s.N), n(s.n), m(s.m) {
  provers_.resize(M);
}

template <class FieldType>
ProverWrapper<FieldType>::~ProverWrapper() {
  for (auto e = 0; e < M; ++e) {
    if (provers_[e])
      delete provers_[e];
  }
}

template <class FieldType>
void ProverWrapper<FieldType>::r1(block &h_gamma) {
  master_seed_.resize(M);

  osuCrypto::PRNG prng;

  // 1
  for (auto e = 0; e < M; ++e) {
    // 1.a
    master_seed_[e].b = osuCrypto::sysRandomSeed();
  }

  // 2 - ** can be parallelized **
//  const size_t nthreads = std::thread::hardware_concurrency();
//  std::vector<std::thread> threads(nthreads);
//
//  for(auto t = 0u; t < nthreads; t++) {
//    threads[t] = std::thread(std::bind(
//        [&](const int bi, const int ei, const int t) {
//          for (auto e = bi; e < ei; ++e) {
//            provers_[e] = new Prover(set_, a_, t_, secret_);
//
//            provers_[e]->r1(master_seed_[e]);
//          }
//        }, t * M / nthreads, (t + 1) == nthreads ? M : (t + 1) * M / nthreads, t));
//  }
//  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });

  for (auto e = 0; e < M; ++e) {
    provers_[e] = new Prover<FieldType>(set_, a_, t_, secret_);

    provers_[e]->r1(master_seed_[e]);
  }

  osuCrypto::SHA1 sha_h_gamma(sizeof(block));
  for (auto e = 0; e < M; ++e) {
    sha_h_gamma.Update(provers_[e]->h_);
  }
  sha_h_gamma.Final(h_gamma_.bytes);

  h_gamma = h_gamma_; // Set out variable
}

template <class FieldType>
void ProverWrapper<FieldType>::r3(const std::vector<bool> &E, std::vector<block> &seed, std::vector<block> &omegaN, block &h_pi) {
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
//  const size_t nthreads = std::thread::hardware_concurrency();
//  std::vector<std::thread> threads(nthreads);

//  for(auto t = 0u; t < nthreads; t++) {
//    threads[t] = std::thread(std::bind(
//        [&](const int bi, const int ei, const int t) {
//          for (auto e = bi; e < ei; ++e) {
//
//          }
//        }, t * M / nthreads, (t + 1) == nthreads ? M : (t + 1) * M / nthreads, t));
//  }
//  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });

  for (auto e = 0; e < M; ++e) {
    if (!E[e]) {
      provers_[e]->r3();
//      std::cout << "Gamma0 " << e << " " << provers_[e]->gamma_[0].halves[0] << " " << provers_[e]->gamma_[0].halves[1] << std::endl;
    }
  }

  // 2.h
  seed.resize(tau);
  omegaN.resize(M - tau);

  osuCrypto::SHA1 sha_h_pi(sizeof(block)); // For step 3
  for (auto e = 0, o_id = 0, s_id = 0; e < M; ++e) {
    if (E[e]) {
      seed[s_id++] = master_seed_[e]; // For out variable
      continue;
    }

    sha_h_pi.Update(provers_[e]->pi_);
    omegaN[o_id++] = provers_[e]->omegaN_; // Set out variable
  }
  sha_h_pi.Final(h_pi_);

  h_pi = h_pi_; // Set out variable
}

template <class FieldType>
void ProverWrapper<FieldType>::r5(const block &seed_ell, block &h_psi) {
  seed_ell_ = seed_ell;

  prng_seed_ell_.SetSeed(seed_ell_.b);

  for (auto e = 0; e < M; ++e) {
    if (E_[e]) {
      continue;
    }

    osuCrypto::SHA1 sha_psi(sizeof(block));

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
  for (auto e = 0; e < M; ++e) {
    if (!E_[e]) {
      provers_[e]->r5();

//      std::cout << e << std::endl;
//      for (auto i = 0; i < N; ++i)
//        std::cout << "\t" << provers_[e]->o_[i] << std::endl;
    }
  }

  osuCrypto::SHA1 sha_h_psi(sizeof(block));

  for (auto e = 0; e < M; ++e) {
    if (!E_[e]) {
      sha_h_psi.Update(provers_[e]->psi_);

      //std::cout << "WE " << e << " " << w_[e_it].halves[0] << " " << w_[e_it].halves[1] << std::endl;
      //std::cout << "PSI " << e << " " << psi_[e_it].halves[0] << psi_[e_it].halves[1] << std::endl;
    }
  }

  sha_h_psi.Final(h_psi_);

  h_psi = h_psi_;
}

template <class FieldType>
void ProverWrapper<FieldType>::r7(const std::vector<int> &i_bar, block &seed_e_bar, std::vector<std::vector<block>> &seed_tree,
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
  for (auto e = 0; e < M; ++e) {
    if (!E_[e])
      provers_[e]->r7(seed_tree[map[e]], gamma_i_bar[map[e]], alpha_i_bar[map[e]],
                      o_i_bar[map[e]], b_square[map[e]], s[map[e]]);
  }
}


}


#endif //LZKP_PROVER_WRAPPER_H
