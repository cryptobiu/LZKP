//
// Created by roee on 9/3/18.
//

#ifndef LZKP_VERIFIER_WRAPPER_H
#define LZKP_VERIFIER_WRAPPER_H


#include "settings.h"
#include "verifier.h"
#include "Mersenne.h"


namespace lzkp {


template <class FieldType>
class VerifierWrapper {
public:
  VerifierWrapper(const Settings &s, const std::vector<std::vector<FieldType>> &a, const std::vector<FieldType> &t);
  ~VerifierWrapper();

  void r2(const block &h_gamma, std::vector<bool> &E);
  void r4(const std::vector<block> &seed, const std::vector<block> &omegaN, const block &h_pi, block &seed_ell);
  void r6(const block &h_psi, std::vector<int> &i_bar);
  bool r8(const block &seed_e_bar, const std::vector<std::vector<block>> &seed_tree, const std::vector<block> &gamma_i_bar,
          const std::vector<std::vector<FieldType>> &alpha_i_bar, const std::vector<FieldType> &o_i_bar, const std::vector<std::vector<FieldType>> &b_square,
          const std::vector<std::vector<FieldType>> &s);

private:
public:
  // Public known values
  const std::vector<std::vector<FieldType>> &a_;
  const std::vector<FieldType> &t_;

  const Settings &set_;
  const int M;
  const int tau;
  const int N;
//  const uint64_t q;
  const int n;
  const int m;

  std::vector<Verifier<FieldType> *> verifiers_;

  block h_gamma_;
  std::vector<bool> E_;

  block h_pi_;
  block seed_ell_;
  osuCrypto::PRNG prng_seed_ell_;

  block h_psi_;

  block seed_e_bar_;
  osuCrypto::PRNG prng_e_bar_;
};


template <class FieldType>
VerifierWrapper<FieldType>::VerifierWrapper(const Settings &s, const std::vector<std::vector<FieldType>> &a, const std::vector<FieldType> &t)
    : a_(a), t_(t), set_(s), M(s.M), tau(s.tau), N(s.N), n(s.n), m(s.m) {
  verifiers_.resize(M);
}

template <class FieldType>
VerifierWrapper<FieldType>::~VerifierWrapper() {
  for (auto e = 0; e < M; ++e) {
    if (verifiers_[e])
      delete verifiers_[e];
  }
}

template <class FieldType>
void VerifierWrapper<FieldType>::r2(const block &h_gamma, std::vector<bool> &E) {
  h_gamma_ = h_gamma;

  E_.resize(M);

  osuCrypto::PRNG prng;
  prng.SetSeed(osuCrypto::sysRandomSeed());

  for (auto i = M - tau + 1; i <= M; ++i) {
    int r = prng.get<osuCrypto::u32>() % i;

    if (!E_[r])
      E_[r] = true;
    else
      E_[i - 1] = true;
  }

  E = E_;
}

template <class FieldType>
void VerifierWrapper<FieldType>::r4(const std::vector<block> &seed, const std::vector<block> &omegaN, const block &h_pi, block &seed_ell) {
//  seed_ = seed;
//  omegaN_ = omegaN;
  h_pi_ = h_pi;

  // 1
  for (auto e = 0, o_id = 0, s_id = 0; e < M; ++e) {
    verifiers_[e] = new Verifier<FieldType>(set_, a_, t_);

    if (E_[e]) {
      // 2.a
      verifiers_[e]->seed_ = seed[s_id++]; //seed_[e_id++];
    }
    else {
      verifiers_[e]->omegaN_ = omegaN[o_id++];
//      std::cout << "e " << e << " " << verifiers_[e]->omegaN_.b << std::endl;
    }
  }

  // ** can be parallelized **
  for (auto e = 0; e < M; ++e) {
    if (E_[e]) {
      verifiers_[e]->r4();
    }
  }

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
void VerifierWrapper<FieldType>::r6(const block &h_psi, std::vector<int> &i_bar) {
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
bool VerifierWrapper<FieldType>::r8(const block &seed_e_bar, const std::vector<std::vector<block>> &seed_tree, const std::vector<block> &gamma_i_bar, const std::vector<std::vector<FieldType>> &alpha_i_bar,
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
  for (auto e = 0; e < M; ++e) {
    if (E_[e])
      continue;

    const int cur_map = map[e];

    verifiers_[e]->r8(seed_tree[cur_map], gamma_i_bar[cur_map], alpha_i_bar[cur_map], o_i_bar[cur_map],
                      b_square[cur_map], s[cur_map]);
  }

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
    //std::cout << "V H " << e << " " << h_[e].halves[0] << " " << h_[e].halves[1] << std::endl;
    sha_h_gamma.Update(verifiers_[e]->h_); // For step 2

    if (E_[e])
      continue;

    sha_h_pi.Update(verifiers_[e]->pi_); // For step 3
    sha_h_psi.Update(verifiers_[e]->psi_); // For step 4
  }
  block h_gamma_computed;
  sha_h_gamma.Final(h_gamma_computed);

  //std::cout << "V h_gammaC " << h_gamma_computed.halves[0] << " " << h_gamma_computed.halves[1] << std::endl;
  //std::cout << "V h_gamma " << h_gamma_.halves[0] << " " << h_gamma_.halves[1] << std::endl;
  if (!eq(h_gamma_.b, h_gamma_computed.b))
    return false;

  // 3
  block pi;
  sha_h_pi.Final(pi);
  //std::cout << "V pi " << pi.halves[0] << " " << pi.halves[1] << std::endl;
  //std::cout << "V h_pi " << h_pi_.halves[0] << " " << h_pi_.halves[1] << std::endl;
  if (!eq(pi.b, h_pi_.b))
    return false;

  // 4
  block psi;
  sha_h_psi.Final(psi);
  //std::cout << "V psi " << psi.halves[0] << " " << psi.halves[1] << std::endl;
  //std::cout << "V h_psi_ " << h_psi_.halves[0] << " " << h_psi_.halves[1] << std::endl;
  if (!eq(psi.b, h_psi_.b))
    return false;

  return true;
}


}


#endif //LZKP_VERIFIER_WRAPPER_H
