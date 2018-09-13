//
// Created by roee on 9/3/18.
//

#ifndef __LZKP_CAC_VERIFIER_LOGIC_H_FILE__
#define __LZKP_CAC_VERIFIER_LOGIC_H_FILE__


#include "parameters.h"
#include "cac_verifier.h"

#include <thread>


namespace lzkp {


template <class FieldType>
class CacVerifierLogic {
public:
  CacVerifierLogic(const Parameters &s, FieldType **&a, FieldType *&t, bool multi_threaded = false);
  ~CacVerifierLogic();

  void r2(const block &h_gamma, std::vector<uint8_t> &E);
  void r4(const std::vector<block> &seed, const std::vector<block> &omegaN, const block &h_pi, block &seed_ell);
  void r6(const block &h_psi, std::vector<int> &i_bar);
  bool r8(const block &seed_e_bar, const std::vector<std::vector<block>> &seed_tree, const std::vector<block> &gamma_i_bar,
          const std::vector<std::vector<FieldType>> &alpha_i_bar, const std::vector<FieldType> &o_i_bar, const std::vector<std::vector<FieldType>> &b_square,
          const std::vector<std::vector<FieldType>> &s);

private:
public:
  // Public known values
  FieldType **&a_;
  FieldType *&t_;

  const Parameters &par_;
  const int M;
  const int tau;
  const int N;
  const int n;
  const int m;

  size_t nthreads_;

  std::vector<CacVerifier<FieldType> *> verifiers_;

  block h_gamma_;
  std::vector<uint8_t> E_;

  block h_pi_;
  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;

  block h_psi_;

  block seed_e_bar_;
  osuCrypto::PRNG prng_e_bar_;
};


template <class FieldType>
CacVerifierLogic<FieldType>::CacVerifierLogic(const Parameters &s, FieldType **&a, FieldType *&t, bool multi_threaded)
    : a_(a), t_(t), par_(s), M(s.M), tau(s.tau), N(s.N), n(s.n), m(s.m) {
  if (multi_threaded)
    nthreads_ = std::thread::hardware_concurrency();
  else
    nthreads_ = 1;

  verifiers_.resize(M);
}

template <class FieldType>
CacVerifierLogic<FieldType>::~CacVerifierLogic() {
  for (auto e = 0; e < M; ++e) {
    if (verifiers_[e])
      delete verifiers_[e];
  }
}

template <class FieldType>
void CacVerifierLogic<FieldType>::r2(const block &h_gamma, std::vector<uint8_t> &E) {
  h_gamma_ = h_gamma;

  E_.resize(M);

  osuCrypto::PRNG prng;
  prng.SetSeed(osuCrypto::sysRandomSeed());

  for (auto i = M - tau + 1; i <= M; ++i) {
    int r = prng.get<osuCrypto::u32>() % i;

    if (!E_[r])
      E_[r] = 0xFF;
    else
      E_[i - 1] = 0xFF;
  }

  E = E_;
}

template <class FieldType>
void CacVerifierLogic<FieldType>::r4(const std::vector<block> &seed, const std::vector<block> &omegaN, const block &h_pi, block &seed_ell) {
  h_pi_ = h_pi;

  // 1
  for (auto e = 0, o_id = 0, s_id = 0; e < M; ++e) {
    verifiers_[e] = new CacVerifier<FieldType>(par_, a_, t_);

    if (E_[e]) {
      // 2.a
      verifiers_[e]->seed_ = seed[s_id++]; //seed_[e_id++];
    }
    else {
      verifiers_[e]->omegaN_ = omegaN[o_id++];
    }
  }

  // ** can be parallelized **
  std::vector<std::thread> threads(nthreads_);

  for(auto t = 0u; t < nthreads_; t++) {
    threads[t] = std::thread(std::bind(
        [&](const int bi, const int ei, const int t) {
          for (auto e = bi; e < ei; ++e) {
            if (E_[e]) {
              verifiers_[e]->r4();
            }
          }
        }, t * M / nthreads_, (t + 1) == nthreads_ ? M : (t + 1) * M / nthreads_, t));
  }
  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });

  // 2
  seed_ell_.b = osuCrypto::sysRandomSeed();
  prng_seed_ell_.SetSeed(seed_ell_.b);

  for (auto e = 0; e < M; ++e) {
    if (!E_[e]) {
      verifiers_[e]->coefficients_.resize(n + m);

      for (auto i = 0; i < n; ++i) {
        verifiers_[e]->coefficients_[i] = FieldType(prng_seed_ell_.get<block>().halves[0]);
      }
      for (auto i = 0; i < m; ++i) {
        verifiers_[e]->coefficients_[n + i] = FieldType(prng_seed_ell_.get<block>().halves[0]);
      }
    }
  }

  seed_ell = seed_ell_;
}

template <class FieldType>
void CacVerifierLogic<FieldType>::r6(const block &h_psi, std::vector<int> &i_bar) {
  h_psi_ = h_psi;

  i_bar.resize(M - tau);

  osuCrypto::PRNG prng;
  prng.SetSeed(osuCrypto::sysRandomSeed());

  for (auto e = 0, e_id = 0; e < M; ++e) {
    if (E_[e])
      continue;

    verifiers_[e]->i_bar_ = prng.get<osuCrypto::u32>() % N;

    i_bar[e_id++] = verifiers_[e]->i_bar_; // Set out variable
  }
}

template <class FieldType>
bool CacVerifierLogic<FieldType>::r8(const block &seed_e_bar, const std::vector<std::vector<block>> &seed_tree, const std::vector<block> &gamma_i_bar, const std::vector<std::vector<FieldType>> &alpha_i_bar,
                         const std::vector<FieldType> &o_i_bar, const std::vector<std::vector<FieldType>> &b_square,
                         const std::vector<std::vector<FieldType>> &s) {
  seed_e_bar_ = seed_e_bar;
  prng_e_bar_.SetSeed(seed_e_bar_.b);

  std::vector<int> map;

  map.resize(M);
  for (auto e = 0, e_it = 0; e < M; ++e) {
    if (!E_[e]) {
      map[e] = e_it;

      e_it++;
    }
  }

  // Generate g_[e] = for step 1.h
  for (auto e = 0; e < M; ++e) {
    if (E_[e])
      continue;

    verifiers_[e]->g_ = prng_e_bar_.get<block>();
  }

  // Generate w_[e] = for step 1.k
  for (auto e = 0; e < M; ++e) {
    if (E_[e])
      continue;

    verifiers_[e]->w_ = prng_e_bar_.get<block>();
  }

  // ** can be parallelized **
  std::vector<std::thread> threads(nthreads_);

  for(auto t = 0u; t < nthreads_; t++) {
    threads[t] = std::thread(std::bind(
        [&](const int bi, const int ei, const int t) {
          for (auto e = bi; e < ei; ++e) {
            if (E_[e])
              continue;

            const int cur_map = map[e];

            verifiers_[e]->r8(seed_tree[cur_map], gamma_i_bar[cur_map], alpha_i_bar[cur_map], o_i_bar[cur_map],
                              b_square[cur_map], s[cur_map]);
          }
        }, t * M / nthreads_, (t + 1) == nthreads_ ? M : (t + 1) * M / nthreads_, t));
  }
  std::for_each(threads.begin(), threads.end(), [](std::thread& x) { x.join(); });

  for (auto e = 0; e < M; ++e) {
    if (E_[e])
      continue;

    if (verifiers_[e]->reject_)
      return false;
  }

  // 2
  osuCrypto::SHA1 sha_h_gamma(sizeof(block));
  osuCrypto::SHA1 sha_h_pi(sizeof(block));
  osuCrypto::SHA1 sha_h_psi(sizeof(block));

  for (auto e = 0; e < M; ++e) {
    sha_h_gamma.Update(verifiers_[e]->h_); // For step 2

    if (E_[e])
      continue;

    sha_h_pi.Update(verifiers_[e]->pi_); // For step 3
    sha_h_psi.Update(verifiers_[e]->psi_); // For step 4
  }

  block h_gamma_computed;
  sha_h_gamma.Final(h_gamma_computed);

  if (!eq(h_gamma_.b, h_gamma_computed.b))
    return false;

  // 3
  block pi;
  sha_h_pi.Final(pi);

  if (!eq(pi.b, h_pi_.b))
    return false;

  // 4
  block psi;
  sha_h_psi.Final(psi);

  if (!eq(psi.b, h_psi_.b))
    return false;

  return true;
}


}


#endif // __LZKP_CAC_VERIFIER_LOGIC_H_FILE__
