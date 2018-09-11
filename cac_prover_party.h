//
// Created by lzkp on 9/7/18.
//

#ifndef __LZKP_CAC_PROVER_PARTY_H_FILE__
#define __LZKP_CAC_PROVER_PARTY_H_FILE__


#include "prover_party.h"

#include "parameters.h"
#include "cac_prover_logic.h"


namespace lzkp {


template <class FieldType>
class CacProverParty : public ProverParty<FieldType> {
public:
  CacProverParty() : ProverParty<FieldType>() {
#ifdef DEBUG
    std::cout << "Constructing CacProverParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl;
#endif
  }
  ~CacProverParty() {
#ifdef DEBUG
    std::cout << "Destructing CacProverParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl;
#endif
  }

  virtual int init(int argc, const char* const argv[]);
  virtual bool runOnline();

  static const int PROTOCOL_TYPE = 0;

protected:
  virtual int parseArguments(int argc, const char* const argv[]);
  virtual int negotiateParameters();
  virtual int generateData();

  // Public known values
  std::vector<std::vector<FieldType>> a_;
  std::vector<FieldType> t_;

  // Prover's secret
  std::vector<FieldType> secret_;

  Parameters par_;
  int M;
  int tau;
  int N;
  int n;
  int m;
  bool is_accepted = false;

  bool multi_threaded_ = false;
};

template<class FieldType>
int CacProverParty<FieldType>::init(int argc, const char* const argv[]) {
  parseArguments(argc, argv);
  int com = this->initCommunication();

  if (com != 0)
    std::cout << "Cannot connect to prover (error " << com << ")" << std::endl;

  int par = this->negotiateParameters();

  if (par != 0)
    std::cout << "Parameters negotiation failed (error " << par << ")" << std::endl;

  return 0;
}

template<class FieldType>
int CacProverParty<FieldType>::parseArguments(int argc, const char* const argv[]) {
  po::variables_map vm;

  po::options_description network("Network options");
  network.add_options()
    ("port", po::value<int>(&this->port_)->required(), "port")
    ;

  po::options_description performence("Performence options");
  performence.add_options()
    ("multi_threaded,x", po::bool_switch(&multi_threaded_), "should execute in multi-threading?")
    ;

  po::options_description parameter("Parameter options");
  parameter.add_options()
    ("M,M", po::value<int>(&M), "number of MPCs 'in the head'")
    ("tau,t", po::value<int>(&tau), "number of games to open")
    ("N,N", po::value<int>(&N), "number of parties in each MPC")
    ("n,n", po::value<int>(&n), "rows in matrix")
    ("m,m", po::value<int>(&m), "columns in matrix")
    ("is_accepted,a", po::bool_switch(&is_accepted), "should proof be accepted?")
    ;

  po::options_description cmdline_options;

  cmdline_options.add(network).add(performence).add(parameter);

  po::store(po::command_line_parser(argc, argv).options(cmdline_options).allow_unregistered().run(), vm);

  notify(vm);

#ifdef DEBUG
  std::cout << "Parsed parameters: " << std::endl;
  std::cout << "\tPort: " << this->port_ << std::endl;
  std::cout << "\tM: " << M << std::endl;
  std::cout << "\ttau: " << tau << std::endl;
  std::cout << "\tN: " << N << std::endl;
  std::cout << "\tn: " << n << std::endl;
  std::cout << "\tm: " << m << std::endl;
  std::cout << "\tis_accepted: " << is_accepted << std::endl;

  if (multi_threaded_)
    std::cout << "\tMulti-threading enabled" << std::endl;
  else
    std::cout << "\tMultit-hreading disabled" << std::endl;
#endif

  return 0;
}

template<class FieldType>
int CacProverParty<FieldType>::negotiateParameters() {
#ifdef DEBUG
  std::cout << "Negotiating protocol parameters..." << std::endl;
#endif

  iovec iov[1];
  ssize_t nwritten, nread;

  int protocol_type;

  iov[0].iov_base = &protocol_type;
  iov[0].iov_len = sizeof(protocol_type);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

#ifdef DEBUG
  std::cout << "\tValidating protocol type... ";
#endif
  assert (protocol_type == CacProverParty::PROTOCOL_TYPE);
#ifdef DEBUG
  std::cout << "done" << std::endl;
#endif

  uint64_t q;

  iov[0].iov_base = &q;
  iov[0].iov_len = sizeof(q);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

  #ifdef DEBUG
  std::cout << "\tValidating field... ";
#endif
  assert (q == (uint64_t)FieldType::p);

  par_ = Parameters(M, tau, N, n, m);

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\tTransmitting protocol parameters... ";
#endif

  iov[0].iov_base = &par_;
  iov[0].iov_len = sizeof(par_);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

#ifdef DEBUG
  std::cout << "done" << std::endl;
#endif

  this->generateData(); // Generate protocol data

#ifdef DEBUG
  std::cout << "\tTransmitting protocol public data... ";
#endif

  iovec *iov2 = new iovec[n + 1];

  for (auto i = 0; i < n; ++i) {
    iov2[i].iov_base = a_[i].data();
    iov2[i].iov_len = a_[i].size() * sizeof(a_[i][0]);
  }
  iov2[n].iov_base = t_.data();
  iov2[n].iov_len = t_.size() * sizeof(t_[0]);
  nwritten = writev(this->sock_, iov2, n + 1);
  assert (nwritten == (int)iov2[0].iov_len * n + (int)iov2[n].iov_len);

  delete[] iov2;

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "Negotiating protocol parameters... done" << std::endl;
#endif

  return 0;
}

template<class FieldType>
int CacProverParty<FieldType>::generateData() {
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a_.resize(n);
  for (auto i = 0; i < n; ++i)
    a_[i].resize(m);

  t_.resize(n);
  secret_.resize(m);

#ifdef DEBUG
  std::cout << "\tStaring data generation..." << std::endl;
  std::cout << "\t\tGenerating matrix (A)... ";
#endif

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
        a_[nn][mm] = FieldType(prng.get<block>().halves[0]);
    }
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\t\tGenerating secret (s)... ";
#endif

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret_[mm] = FieldType(prng.get<block>().bytes[0] % 2);
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\t\tCalculating vector (t)... ";
#endif

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t_[nn] = FieldType(0);

    for (auto mm = 0; mm < m; ++mm) {
      t_[nn] += a_[nn][mm] * secret_[mm];
    }
  }

  if (!is_accepted) {
#ifdef DEBUG
    std::cout << "done" << std::endl;
    std::cout << "\t\tProof should not be accepted, adding noise to the secret... ";
#endif
    secret_[0] += FieldType(1);
  }

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\tData generation done" << std::endl;
#endif

  return 0;
}

template<class FieldType>
bool CacProverParty<FieldType>::runOnline() {
  CacProverLogic<FieldType> p(par_, a_, t_, secret_, multi_threaded_);

  iovec *iov = new iovec[(M - tau) * 4 + 3]; // 1 + (M - tau) + 1 + (M - tau) + 1 + i_id, i_id maximum value is 2 * (M - tau)
  ssize_t nwritten, nread;

  // ** Round 1 **
  block h_gamma;
  p.r1(h_gamma); // Run round 1
  iov[0].iov_base = &h_gamma;
  iov[0].iov_len = sizeof(h_gamma);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  // ** Round 2 output **
  std::vector<uint8_t> E(M);
  iov[0].iov_base = E.data();
  iov[0].iov_len = E.size() * sizeof(E[0]);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

  // ** Round 3 **
  std::vector<block> seed, omegaN;
  block h_pi;
  p.r3(E, seed, omegaN, h_pi); // Run round 3
  iov[0].iov_base = seed.data();
  iov[0].iov_len = seed.size() * sizeof(seed[0]);
  iov[1].iov_base = omegaN.data();
  iov[1].iov_len = omegaN.size() * sizeof(omegaN[0]);
  iov[2].iov_base = &h_pi;
  iov[2].iov_len = sizeof(h_pi);
  nwritten = writev(this->sock_, iov, 3);
  assert (nwritten == (int)(iov[0].iov_len + iov[1].iov_len + iov[2].iov_len));

  // ** Round 4 output **
  block seed_ell;
  iov[0].iov_base = &seed_ell;
  iov[0].iov_len = sizeof(seed_ell);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

  // ** Round 5 **
  block h_psi;
  p.r5(seed_ell, h_psi); // Run round 5
  iov[0].iov_base = &h_psi;
  iov[0].iov_len = sizeof(h_psi);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  // ** Round 6 output **
  std::vector<int> i_bar(M - tau);
  iov[0].iov_base = i_bar.data();
  iov[0].iov_len = i_bar.size() * sizeof(i_bar[0]);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

  // ** Round 7 **
  block seed_e_bar;
  std::vector<std::vector<block>> seed_tree;
  std::vector<block> gamma_i_bar;
  std::vector<std::vector<FieldType>> alpha_i_bar, b_square, s;
  std::vector<FieldType> o_i_bar;
  p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
  auto iov_id = 0;
  iov[iov_id].iov_base = &seed_e_bar;
  iov[iov_id++].iov_len = sizeof(seed_e_bar);
  for (auto i = 0; i < M - tau; ++i) {
    iov[iov_id].iov_base = seed_tree[i].data();
    iov[iov_id++].iov_len = seed_tree[i].size() * sizeof(seed_tree[i][0]);
  }
  iov[iov_id].iov_base = gamma_i_bar.data();
  iov[iov_id++].iov_len = gamma_i_bar.size() * sizeof(gamma_i_bar[0]);
  for (auto i = 0; i < M - tau; ++i) {
    iov[iov_id].iov_base = alpha_i_bar[i].data();
    iov[iov_id++].iov_len = alpha_i_bar[i].size() * sizeof(alpha_i_bar[i][0]);
  }
  iov[iov_id].iov_base = o_i_bar.data();
  iov[iov_id++].iov_len = o_i_bar.size() * sizeof(o_i_bar[0]);
  int i_id = 0;
  for (auto e = 0, e_id = 0; e < M; ++e) {
    if (E[e])
      continue;

    if (p.provers_[e]->i_bar_ != N - 1) {
      iov[iov_id].iov_base = b_square[e_id].data();
      iov[iov_id++].iov_len = b_square[e_id].size() * sizeof(b_square[e_id][0]);
      iov[iov_id].iov_base = s[e_id].data();
      iov[iov_id++].iov_len = s[e_id].size() * sizeof(s[e_id][0]);

      i_id += 2;
    }

    e_id++;
  }
  nwritten = writev(this->sock_, iov, 1 + (M - tau) + 1 + (M - tau) + 1 + i_id);
  std::cout << nwritten << std::endl;
  std::cout << (int)(iov[0].iov_len + iov[1].iov_len * (M - tau) + iov[M - tau + 1].iov_len +
                     iov[M - tau + 2].iov_len * (M - tau) + iov[2 * (M - tau) + 2].iov_len +
                     iov[2 * (M - tau) + 3].iov_len * i_id) << std::endl;
  assert (nwritten == (int)(iov[0].iov_len + iov[1].iov_len * (M - tau) + iov[M - tau + 1].iov_len +
                            iov[M - tau + 2].iov_len * (M - tau) + iov[2 * (M - tau) + 2].iov_len +
                            iov[2 * (M - tau) + 3].iov_len * i_id));

  delete[] iov;

  return true;
}


}


#endif // __LZKP_CAC_PROVER_PARTY_H_FILE__
