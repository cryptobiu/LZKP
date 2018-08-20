//
// Created by roee on 7/26/18.
//

#include "verifier.h"

#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Crypto/sha1.h>

#include "seedtree.h"

using namespace lzkp;

Verifier::Verifier(const Settings &s) : M(s.M), N(s.N), q(s.q), m(s.m), n(s.n), tau(s.tau) {
  E_.resize(M);

  NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***
}

const std::vector<bool> &Verifier::r2(const block &h_gamma) {
  h_gamma_ = h_gamma;

  osuCrypto::PRNG prng;

  prng.SetSeed(osuCrypto::sysRandomSeed());

  for (auto i = M - tau + 1; i <= M; ++i) {
    int r = prng.get<osuCrypto::u32>() % i;

    if (!E_[r])
      E_[r] = true;
    else
      E_[i - 1] = true;
  }

  return E_;
}

void Verifier::r4(const std::vector<block> &seed, const std::vector<block> &omega, const block &h_pi, std::vector<std::vector<NTL::ZZ_p>> &coefficients) {
  seed_ = seed;
  omega_ = omega;
  h_pi_ = h_pi;

  coefficients.resize(M - tau);

  seed_tree_.resize(M);
  r_.resize(M);
  b_.resize(M);
  b_square_.resize(M);
  gamma_.resize(M);
  h_.resize(M);

  int seed_id = 0;
  int coef_id = 0;

  for (auto e = 0; e < M; ++e) {
    if (E_[e]) { // 1
      // *1* - 1.a
      seed_tree_[e].resize(N);
      gamma_[e].resize(N);

      seed_tree_[e].generate(seed[seed_id++]); // Get the relevant seed

      r_[e].resize(N);
      b_[e].resize(m);
      b_square_[e].resize(m);

      for (auto mm = 0; mm < m; ++mm){
        b_[e][mm].resize(N);
        b_square_[e][mm].resize(N);
      }

      // *1* - 1.b
      for (auto i = 0; i < N - 1; ++i) {
        r_[e][i] = seed_tree_[e].getBlock(i);

        for (auto mm = 0; mm < m; ++mm) {
          b_[e][mm][i]        = NTL::ZZ_p(seed_tree_[e].getBlock(i).halves[0]); // NEED TO FIND A WAY TO USE ALL 128 BITS
          b_square_[e][mm][i] = NTL::ZZ_p(seed_tree_[e].getBlock(i).halves[0]);
        }
      }

      // *1* - 1.c
      r_[e][N-1] = seed_tree_[e].getBlock(N - 1);

      for (auto mm = 0; mm < m; ++mm) {
        b_[e][mm][N-1]        = NTL::ZZ_p(seed_tree_[e].getBlock(N - 1).halves[0]);
      }

      // *1* - 1.d
      for (auto k = 0; k < m; ++k) {
        b_square_[e][k][N - 1] = NTL::ZZ_p(0);
        for (auto i = 0; i < N; ++i) {
          b_square_[e][k][N - 1] += b_[e][k][i];
        }
        b_square_[e][k][N - 1] *= b_square_[e][k][N - 1];
        for (auto i = 0; i < N - 1; ++i) {
          b_square_[e][k][N - 1] -= b_square_[e][k][i];
        }
      }

      // *1* - 1.e
      osuCrypto::SHA1 sha_gamma(sizeof(block));
      for (auto i = 0; i < N - 1; ++i) {
        block blk = seed_tree_[e].getSeed(i);
        sha_gamma.Reset();
        sha_gamma.Update(blk);
        sha_gamma.Update(r_[e][i]);
        sha_gamma.Final(gamma_[e][i]);
      }
      block blk = seed_tree_[e].getSeed(N - 1);
      sha_gamma.Reset();
      sha_gamma.Update(blk);
      unsigned char buf[1024]; // MB NEED TO ZERO BUFFER, OR TO HASH ONLY PART OF IT
      for (auto i = 0; i < m; ++i) {
        NTL::BytesFromZZ(buf, b_square_[e][i][N - 1]._ZZ_p__rep, NTL::NumBytes(b_square_[e][i][N - 1]._ZZ_p__rep));
        sha_gamma.Update(buf, NTL::NumBytes(b_square_[e][i][N - 1]._ZZ_p__rep));
      }
      sha_gamma.Update(r_[e][N - 1]);
      sha_gamma.Final(gamma_[e][N - 1]);

      // *1* - 1.f
      osuCrypto::SHA1 sha_h(sizeof(block));
      for (auto i = 0; i < N; ++i) {
        sha_h.Update(gamma_[e][i]);
      }
      sha_h.Final(h_[e]); // mb need to zero it first
    }
    else { // 2
      osuCrypto::PRNG prng;
      prng.SetSeed(osuCrypto::sysRandomSeed());

      coefficients[coef_id].resize(n + m);

      for (auto i = 0; i < n; ++i) {
        coefficients[coef_id][i] = NTL::ZZ_p(prng.get<block>().halves[0]);
      }
      for (auto i = 0; i < m; ++i) {
        coefficients[coef_id][n + i] = NTL::ZZ_p(prng.get<block>().halves[0]);
      }

      coef_id++;
    }
  }
}

void Verifier::r6(const block &h_psi, std::vector<int> &i_bar) {
  h_psi_ = h_psi;

  i_bar.resize(M - tau);

  osuCrypto::PRNG prng;

  prng.SetSeed(osuCrypto::sysRandomSeed());

  for (auto e = 0; e < M - tau; ++e) {
    i_bar[e] = prng.get<osuCrypto::u32>() % N;
  }
}
