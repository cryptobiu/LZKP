//
// Created by lzkp on 9/7/18.
//

#ifndef __LZKP_SAC_PROVER_PARTY_H_FILE__
#define __LZKP_SAC_PROVER_PARTY_H_FILE__


#include "prover_party.h"

#include "parameters.h"
#include "sac_prover_logic.h"


namespace lzkp {


template <class FieldType>
class SacProverParty : public ProverParty<FieldType> {
public:
  SacProverParty() : ProverParty<FieldType>() {
    debug("Constructing SacProverParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl);
  }
  ~SacProverParty() {
    debug("Destructing SacProverParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl);
  }

  virtual int init(int argc, const char* const argv[]);
  virtual bool runOnline();

  static const int PROTOCOL_TYPE = 1;

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
  int N;
  int n;
  int m;
  bool is_accepted = false;

  bool multi_threaded_ = false;
};

template<class FieldType>
int SacProverParty<FieldType>::init(int argc, const char* const argv[]) {
  parseArguments(argc, argv);
  int com = this->initCommunication();

  if (com != 0)
    std::cerr << "Cannot connect to prover (error " << com << ")" << std::endl;

  int par = this->negotiateParameters();

  if (par != 0)
    std::cerr << "Parameters negotiation failed (error " << par << ")" << std::endl;

  return 0;
}

template<class FieldType>
int SacProverParty<FieldType>::parseArguments(int argc, const char* const argv[]) {
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
    ("N,N", po::value<int>(&N), "number of parties in each MPC")
    ("n,n", po::value<int>(&n), "rows in matrix")
    ("m,m", po::value<int>(&m), "columns in matrix")
    ("is_accepted,a", po::bool_switch(&is_accepted), "should proof be accepted?")
    ;

  po::options_description cmdline_options;

  cmdline_options.add(network).add(performence).add(parameter);

  po::store(po::command_line_parser(argc, argv).options(cmdline_options).allow_unregistered().run(), vm);

  notify(vm);

  debug("Parsed parameters: " << std::endl);
  debug("\tPort: " << this->port_ << std::endl);
  debug("\tM: " << M << std::endl);
  debug("\tN: " << N << std::endl);
  debug("\tn: " << n << std::endl);
  debug("\tm: " << m << std::endl);
  debug("\tis_accepted: " << is_accepted << std::endl);

  if (multi_threaded_)
    debug("\tMulti-threading enabled" << std::endl);
  else
  debug("\tMultit-hreading disabled" << std::endl);

  return 0;
}

template<class FieldType>
int SacProverParty<FieldType>::negotiateParameters() {
  debug("Negotiating protocol parameters..." << std::endl);

  iovec iov[1];
  ssize_t nwritten, nread;

  int protocol_type;

  iov[0].iov_base = &protocol_type;
  iov[0].iov_len = sizeof(protocol_type);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

  debug("\tValidating protocol type... ");

  assert (protocol_type == SacProverParty::PROTOCOL_TYPE);

  debug("done" << std::endl);

  uint64_t q;

  iov[0].iov_base = &q;
  iov[0].iov_len = sizeof(q);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

  debug("\tValidating field... ");

  assert (q == (uint64_t)FieldType::p);

  par_ = Parameters(M, 0, N, n, m);

  debug("done" << std::endl);
  debug("\tTransmitting protocol parameters... ");

  iov[0].iov_base = &par_;
  iov[0].iov_len = sizeof(par_);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  debug("done" << std::endl);

  this->generateData(); // Generate protocol data

  debug("\tTransmitting protocol public data... ");

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

  debug("done" << std::endl);
  debug("Negotiating protocol parameters... done" << std::endl);

  return 0;
}

template<class FieldType>
int SacProverParty<FieldType>::generateData() {
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a_.resize(n);
  for (auto i = 0; i < n; ++i)
    a_[i].resize(m);

  t_.resize(n);
  secret_.resize(m);

  debug("\tStaring data generation..." << std::endl);
  debug("\t\tGenerating matrix (A)... ");

  // Random matrix A
  for (auto nn = 0; nn < n; ++nn) {
    for (auto mm = 0; mm < m; ++mm) {
        a_[nn][mm] = FieldType(prng.get<block>().halves[0]);
    }
  }

  debug("done" << std::endl);
  debug("\t\tGenerating secret (s)... ");

  // Random vector secret
  for (auto mm = 0; mm < m; ++mm) {
    secret_[mm] = FieldType(prng.get<block>().bytes[0] % 2);
  }

  debug("done" << std::endl);
  debug("\t\tCalculating vector (t)... ");

  // Calculate t
  for (auto nn = 0; nn < n; ++nn) {
    t_[nn] = FieldType(0);

    for (auto mm = 0; mm < m; ++mm) {
      t_[nn] += a_[nn][mm] * secret_[mm];
    }
  }

  if (!is_accepted) {
    debug("done" << std::endl);
    debug("\t\tProof should not be accepted, adding noise to the secret... ");

    secret_[0] += FieldType(1);
  }

  debug("done" << std::endl);
  debug("\tData generation done" << std::endl);

  return 0;
}

template<class FieldType>
bool SacProverParty<FieldType>::runOnline() {
  SacProverLogic<FieldType> p(par_, a_, t_, secret_, multi_threaded_);

  iovec *iov = new iovec[M * 5 + 4]; // 1 + M + 1 + M + 2 + i_id, i_id maximum value is 3 * M
  ssize_t nwritten, nread;

  // ** Round 1 **
  block h_gamma;
  p.r1(h_gamma); // Run round 1
  iov[0].iov_base = &h_gamma;
  iov[0].iov_len = sizeof(h_gamma);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  // ** Round 2 output **
  block seed_ell;
  iov[0].iov_base = &seed_ell;
  iov[0].iov_len = sizeof(seed_ell);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

  // ** Round 3 **
  block h_pi, h_psi, h_theta;
  p.r3(seed_ell, h_pi, h_psi, h_theta); // Run round 3
  iov[0].iov_base = &h_pi;
  iov[0].iov_len = sizeof(h_pi);
  iov[1].iov_base = &h_psi;
  iov[1].iov_len = sizeof(h_psi);
  iov[2].iov_base = &h_theta;
  iov[2].iov_len = sizeof(h_theta);
  nwritten = writev(this->sock_, iov, 3);
  assert (nwritten == (int)(iov[0].iov_len + iov[1].iov_len + iov[2].iov_len));

  // ** Round 4 output **
  std::vector<int> i_bar(M);
  iov[0].iov_base = i_bar.data();
  iov[0].iov_len = i_bar.size() * sizeof(i_bar[0]);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

  // ** Round 5 **
  block seed_global;
  std::vector<std::vector<block>> seed_tree;
  std::vector<block> gamma_i_bar;
  std::vector<std::vector<FieldType>> alpha_i_bar, b_square, s, s_square;
  std::vector<FieldType> o_i_bar, v_i_bar;
  p.r5(i_bar, seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 5
  auto iov_id = 0;
  iov[iov_id].iov_base = &seed_global;
  iov[iov_id++].iov_len = sizeof(seed_global);
  for (auto i = 0; i < M; ++i) {
    iov[iov_id].iov_base = seed_tree[i].data();
    iov[iov_id++].iov_len = seed_tree[i].size() * sizeof(seed_tree[i][0]);
  }
  iov[iov_id].iov_base = gamma_i_bar.data();
  iov[iov_id++].iov_len = gamma_i_bar.size() * sizeof(gamma_i_bar[0]);
  for (auto i = 0; i < M; ++i) {
    iov[iov_id].iov_base = alpha_i_bar[i].data();
    iov[iov_id++].iov_len = alpha_i_bar[i].size() * sizeof(alpha_i_bar[i][0]);
  }
  iov[iov_id].iov_base = o_i_bar.data();
  iov[iov_id++].iov_len = o_i_bar.size() * sizeof(o_i_bar[0]);
  iov[iov_id].iov_base = v_i_bar.data();
  iov[iov_id++].iov_len = v_i_bar.size() * sizeof(v_i_bar[0]);
  int i_id = 0;
  for (auto e = 0; e < M; ++e) {
    if (p.provers_[e]->i_bar_ != N - 1) {
      iov[iov_id].iov_base = b_square[e].data();
      iov[iov_id++].iov_len = b_square[e].size() * sizeof(b_square[e][0]);
      iov[iov_id].iov_base = s[e].data();
      iov[iov_id++].iov_len = s[e].size() * sizeof(s[e][0]);
      iov[iov_id].iov_base = s_square[e].data();
      iov[iov_id++].iov_len = s_square[e].size() * sizeof(s_square[e][0]);

      i_id += 3;
    }
  }
  nwritten = writev(this->sock_, iov, 1 + M + 1 + M + 2 + i_id);
  assert (nwritten == (int)(iov[0].iov_len + iov[1].iov_len * M + iov[M + 1].iov_len +
                            iov[M + 2].iov_len * M + iov[2 * M + 2].iov_len + iov[2 * M + 3].iov_len +
                            iov[2 * M + 4].iov_len * i_id));

  delete[] iov;

  return true;
}


}


#endif // __LZKP_SAC_PROVER_PARTY_H_FILE__
