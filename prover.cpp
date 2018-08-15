#include "prover.h"
#include "seedtree.h"

using namespace lzkp;

Prover::Prover(const Settings &s) : M(s.M), N(s.N), q(s.q), m(s.m), t(s.t) {
  st_.resize(M);
  master_seed_.resize(M);
  r_.resize(M);
  b_.resize(M);
  b_square_.resize(M);
  gamma_.resize(M);
}

block Prover::r1() {
  osuCrypto::PRNG prng;
  NTL::ZZ_p::init(NTL::ZZ(q));

  osuCrypto::SHA1 h_gamma(sizeof(block));

  for (auto e = 0; e < M; ++e) {
    st_[e].Resize(N);
    gamma_[e].resize(N);

    // 1.a
    master_seed_[e].b = osuCrypto::sysRandomSeed();
    st_[e].Generate(master_seed_[e]);

    r_[e].resize(N);
    b_[e].resize(m);
    b_square_[e].resize(m);

    for (auto mm = 0; mm < m; ++mm){
      b_[e][mm].resize(N);
      b_square_[e][mm].resize(N);
    }

    for (auto i = 0; i < N - 1; ++i) {
      block blk = st_[e][N - 1 + i];
      prng.SetSeed(blk.b);

      // 1.b
      r_[e][i] = prng.get<block>();

      for (auto mm = 0; mm < m; ++mm) {
        b_[e][mm][i]        = NTL::ZZ(prng.get<block>().halves[0]);
        b_square_[e][mm][i] = NTL::ZZ(prng.get<block>().halves[0]);
      }
    }

    // 1.c
    block blk = st_[e][2 * N - 1];
    r_[e][N-1] = prng.get<block>();

    for (auto mm = 0; mm < m; ++mm) {
      b_[e][mm][N-1]        = NTL::ZZ(prng.get<block>().halves[0]);
    }

    // 1.d
    for (auto k = 0; k < m; ++k) {
      b_square_[e][k][N - 1] = NTL::ZZ(0);
      for (auto i = 0; i < N; ++i) {
        b_square_[e][k][N - 1] += b_[e][k][i];
      }
      b_square_[e][k][N - 1] *= b_square_[e][k][N - 1];
      for (auto i = 0; i < N - 1; ++i) {
        b_square_[e][k][N - 1] -= b_square_[e][k][i];
      }
    }

    // 1.e
    osuCrypto::SHA1 sha(sizeof(block));
    for (auto i = 0; i < N - 1; ++i) {
      blk = st_[e][N - 1 + i];
      sha.Reset();
      sha.Update(&blk, 1);
      sha.Update(&r_[e][i], 1);
      sha.Final(gamma_[e][i]);
    }
    blk = st_[e][2 * N - 2];
    sha.Reset();
    sha.Update(&blk, 1);
    unsigned char buf[200];
    for (auto i = 0; i < m; ++i) {
      NTL::BytesFromZZ(buf, b_square_[e][i][N - 1], NTL::NumBytes(b_square_[e][i][N - 1]));
      sha.Update(buf, NTL::NumBytes(b_square_[e][i][N - 1]));
    }
    sha.Update(&r_[e][N - 1], 1);
    sha.Final(gamma_[e][N - 1]);

    // 1.f
    osuCrypto::SHA1 h_e(sizeof(block));
    for (auto i = 0; i < N - 1; ++i) {
      h_e.Update(&gamma_[e][i], 1);
    }
    h_e.Final(blk.bytes); // mb need to zero it first

    h_gamma.Update(&blk, 1);
  }

  block blk;
  h_gamma.Final(blk.bytes);

  return blk;
}

Prover::~Prover() {
//  for (auto e = 0; e < M; ++e) {
//    for (auto i = 0; i < N; ++i) {
//      delete gamma_[e][i];
//    }
//  }
}

void Prover::r3(std::vector<bool> E) {
  // 1
  auto seed_e_bar = osuCrypto::sysRandomSeed();

  osuCrypto::PRNG prng;
  prng.SetSeed(seed_e_bar);

  // 2
  for (auto e = 0; e < M; ++e) {
    if (E[e])
      continue;

    // 2.a
    auto g_e = prng.get<block>();

    // 2.b
    for (auto i = 0; i < N - 1; ++i) {

    }

  }
}
