//
// Created by lzkp on 9/7/18.
//

#ifndef __LZKP_SAC_VERIFIER_PARTY_H_FILE__
#define __LZKP_SAC_VERIFIER_PARTY_H_FILE__


#include "verifier_party.h"

#include "parameters.h"
#include "sac_verifier_logic.h"


namespace lzkp {


template <class FieldType>
class SacVerifierParty : public VerifierParty<FieldType> {
public:
  SacVerifierParty() : VerifierParty<FieldType>() {
    debug("Constructing SacVerifierParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl);
  }
  ~SacVerifierParty() {
    debug("Destructing SacVerifierParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl);
  }

  virtual int init(int argc, const char* const argv[]);
  virtual bool runOnline();

  static const int PROTOCOL_TYPE = 1;

protected:
  virtual int parseArguments(int argc, const char* const argv[]);
  virtual int negotiateParameters();

  // Public known values
  std::vector<std::vector<FieldType>> a_;
  std::vector<FieldType> t_;

  Parameters par_;

  bool multi_threaded_ = false;

};

template<class FieldType>
int SacVerifierParty<FieldType>::init(int argc, const char* const argv[]) {
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
int SacVerifierParty<FieldType>::parseArguments(int argc, const char* const argv[]) {
  po::variables_map vm;

  po::options_description network("Network options");
  network.add_options()
    ("ip", po::value<std::string>(&this->ip_)->required(), "other party IP (required only for verifier)")
    ("port", po::value<int>(&this->port_)->required(), "port")
    ;

  po::options_description performence("Performence options");
  performence.add_options()
    ("multi_threaded,x", po::bool_switch(&multi_threaded_), "should execute in multi-threading?")
    ;

  po::options_description cmdline_options;

  cmdline_options.add(network).add(performence);

  po::store(po::command_line_parser(argc, argv).options(cmdline_options).allow_unregistered().run(), vm);

  notify(vm);

  debug("Parsed parameters: " << std::endl);
  debug("\tIP: " << this->ip_ << std::endl);
  debug("\tPort: " << this->port_ << std::endl);

  if (multi_threaded_)
    debug("\tMulti-threading enabled" << std::endl);
  else
    debug("\tMulti-threading disabled" << std::endl);

  return 0;
}

template<class FieldType>
int SacVerifierParty<FieldType>::negotiateParameters() {
  debug("Negotiating protocol parameters..." << std::endl);

  iovec iov[1];
  ssize_t nwritten, nread;

  int protocol_type = SacVerifierParty::PROTOCOL_TYPE;

  debug("\tTransmitting protocol type... ");

  iov[0].iov_base = &protocol_type;
  iov[0].iov_len = sizeof(protocol_type);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  debug("done" << std::endl);
  debug("\tTransmitting field characteristic ... ");


  uint64_t q = FieldType::p;

  iov[0].iov_base = &q;
  iov[0].iov_len = sizeof(q);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  debug("done" << std::endl);
  debug("\tReceiving protocol parameters... ");

  iov[0].iov_base = &this->par_;
  iov[0].iov_len = sizeof(this->par_);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

  debug("done" << std::endl);
  debug("\t\tM: " << this->par_.M << std::endl);
  debug("\t\tN: " << this->par_.N << std::endl);
  debug("\t\tn: " << this->par_.n << std::endl);
  debug("\t\tm: " << this->par_.m << std::endl);
  debug("\tReceiving protocol public data... ");

  iovec *iov2 = new iovec[par_.n + 1];

  a_.resize(par_.n);
  for (auto i = 0; i < par_.n; ++i) {
    a_[i].resize(par_.m);
    iov2[i].iov_base = a_[i].data();
    iov2[i].iov_len = a_[i].size() * sizeof(a_[i][0]);
  }
  t_.resize(par_.n);
  iov2[par_.n].iov_base = t_.data();
  iov2[par_.n].iov_len = t_.size() * sizeof(t_[0]);
  nread = readv(this->sock_, iov2, par_.n + 1);
  assert (nread == (int)iov2[0].iov_len * par_.n + (int)iov2[par_.n].iov_len);

  delete[] iov2;

  debug("done" << std::endl);
  debug("Negotiating protocol parameters... done" << std::endl);

  return 0;
}

template<class FieldType>
bool SacVerifierParty<FieldType>::runOnline() {
  SacVerifierLogic<FieldType> v(par_, a_, t_, multi_threaded_);

  iovec *iov = new iovec[par_.M * 5 + 4]; // 1 + par_.M + 1 + par_.M + 2 + par_.M * 3
  ssize_t nwritten, nread;

  // ** Round 1 output **
  block h_gamma;
  iov[0].iov_base = &h_gamma;
  iov[0].iov_len = sizeof(h_gamma);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

  // ** Round 2 **
  block seed_ell;
  v.r2(h_gamma, seed_ell); // Run round 2
  iov[0].iov_base = &seed_ell;
  iov[0].iov_len = sizeof(seed_ell);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  // ** Round 3 output **
  block h_pi, h_psi, h_theta;
  iov[0].iov_base = &h_pi;
  iov[0].iov_len = sizeof(h_pi);
  iov[1].iov_base = &h_psi;
  iov[1].iov_len = sizeof(h_psi);
  iov[2].iov_base = &h_theta;
  iov[2].iov_len = sizeof(h_theta);
  nread = readv(this->sock_, iov, 3);
  assert (nread == (int)(iov[0].iov_len + iov[1].iov_len + iov[2].iov_len));

  // ** Round 4 **
  std::vector<int> i_bar;
  v.r4(h_pi, h_psi, h_theta, i_bar); // Run round 4
  iov[0].iov_base = i_bar.data();
  iov[0].iov_len = i_bar.size() * sizeof(i_bar[0]);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  // ** Round 5 output **
  block seed_global;
  std::vector<std::vector<block>> seed_tree(par_.M);
  std::vector<block> gamma_i_bar(par_.M);
  std::vector<std::vector<FieldType>> alpha_i_bar(par_.M), b_square(par_.M), s(par_.M), s_square(par_.M);
  std::vector<FieldType> o_i_bar(par_.M), v_i_bar(par_.M);
  auto iov_id = 0;
  iov[iov_id].iov_base = &seed_global;
  iov[iov_id++].iov_len = sizeof(seed_global);
  for (auto i = 0; i < par_.M; ++i) {
    seed_tree[i].resize(log(par_.N)/log(2));
    iov[iov_id].iov_base = seed_tree[i].data();
    iov[iov_id++].iov_len = seed_tree[i].size() * sizeof(seed_tree[i][0]);
  }
  iov[iov_id].iov_base = gamma_i_bar.data();
  iov[iov_id++].iov_len = gamma_i_bar.size() * sizeof(gamma_i_bar[0]);
  for (auto i = 0; i < par_.M; ++i) {
    alpha_i_bar[i].resize(par_.m);
    iov[iov_id].iov_base = alpha_i_bar[i].data();
    iov[iov_id++].iov_len = alpha_i_bar[i].size() * sizeof(alpha_i_bar[i][0]);
  }
  iov[iov_id].iov_base = o_i_bar.data();
  iov[iov_id++].iov_len = o_i_bar.size() * sizeof(o_i_bar[0]);
  iov[iov_id].iov_base = v_i_bar.data();
  iov[iov_id++].iov_len = v_i_bar.size() * sizeof(v_i_bar[0]);
  int i_id = 0;
  for (auto e = 0; e < par_.M; ++e) {
    if (v.verifiers_[e]->i_bar_ != par_.N - 1) {
      b_square[e].resize(par_.m);
      iov[iov_id].iov_base = b_square[e].data();
      iov[iov_id++].iov_len = b_square[e].size() * sizeof(b_square[e][0]);
      s[e].resize(par_.m);
      iov[iov_id].iov_base = s[e].data();
      iov[iov_id++].iov_len = s[e].size() * sizeof(s[e][0]);
      s_square[e].resize(par_.m);
      iov[iov_id].iov_base = s_square[e].data();
      iov[iov_id++].iov_len = s_square[e].size() * sizeof(s_square[e][0]);

      i_id += 3;
    }
  }
  nread = readv(this->sock_, iov, 1 + par_.M + 1 + par_.M + 2 + i_id);
  assert (nread == (int)(iov[0].iov_len + iov[1].iov_len * par_.M + iov[par_.M + 1].iov_len +
                            iov[par_.M + 2].iov_len * par_.M + iov[2 * par_.M + 2].iov_len + iov[2 * par_.M + 3].iov_len +
                            iov[2 * par_.M + 4].iov_len * i_id));


  bool flag = v.r6(seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 6

  if (flag)
    std::cout << "Proof accepted" << std::endl;
  else
    std::cout << "Proof rejected" << std::endl;

  delete[] iov;

  return true;
}


}


#endif // __LZKP_SAC_VERIFIER_PARTY_H_FILE__