//
// Created by lzkp on 9/7/18.
//

#ifndef __LZKP_CAC_VERIFIER_PARTY_H_FILE__
#define __LZKP_CAC_VERIFIER_PARTY_H_FILE__


#include "verifier_party.h"

#include "settings.h"
#include "cac_verifier_logic.h"


namespace lzkp {


template <class FieldType>
class CacVerifierParty : public VerifierParty<FieldType> {
public:
  CacVerifierParty() : VerifierParty<FieldType>() {
#ifdef DEBUG
    std::cout << "Constructing CacVerifierParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl;
#endif
  }
  ~CacVerifierParty() {
#ifdef DEBUG
    std::cout << "Destructing CacVerifierParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl;
#endif
  }

  virtual int init(int argc, const char* const argv[]);
  virtual bool runOnline();

  static const int PROTOCOL_TYPE = 0;

protected:
  virtual int parseArguments(int argc, const char* const argv[]);
  virtual int negotiateParameters();

  // Public known values
  std::vector<std::vector<FieldType>> a_;
  std::vector<FieldType> t_;

  Settings set_;
};

template<class FieldType>
int CacVerifierParty<FieldType>::init(int argc, const char* const argv[]) {
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
int CacVerifierParty<FieldType>::parseArguments(int argc, const char* const argv[]) {
  po::variables_map vm;

  po::options_description network("Network options");
  network.add_options()
    ("ip", po::value<std::string>(&this->ip_)->required(), "other party IP (required only for verifier)")
    ("port", po::value<int>(&this->port_)->required(), "port")
    ;

//  po::options_description parameter("Parameter options");
//  parameter.add_options()
//          ("q,q", po::value<int>(&q), "field type (15=15BitField, 31=MersenneInt, 61=MersenneLong)")
//          ("M,M", po::value<int>(&M), "number of MPCs 'in the head'")
//          ("tau,t", po::value<int>(&tau), "number of games to open")
//          ("N,N", po::value<int>(&N), "number of parties in each MPC")
//          ("n,n", po::value<int>(&n), "rows in matrix")
//          ("m,m", po::value<int>(&m), "columns in matrix")
//          ("is_accepted,a", po::bool_switch(&is_accepted), "should proof be accepted?")
//          ;

  po::options_description cmdline_options;

  cmdline_options.add(network);

  po::store(po::command_line_parser(argc, argv).options(cmdline_options).allow_unregistered().run(), vm);

  notify(vm);

#ifdef DEBUG
  std::cout << "Parsed parameters: " << std::endl;
  std::cout << "\tIP: " << this->ip_ << std::endl;
  std::cout << "\tPort: " << this->port_ << std::endl;
#endif

  return 0;
}

template<class FieldType>
int CacVerifierParty<FieldType>::negotiateParameters() {
#ifdef DEBUG
  std::cout << "Negotiating protocol parameters..." << std::endl;
#endif

  iovec iov[1];
  ssize_t nwritten, nread;

  int protocol_type = CacVerifierParty::PROTOCOL_TYPE;

#ifdef DEBUG
  std::cout << "\tTransmitting protocol type... ";
#endif
  iov[0].iov_base = &protocol_type;
  iov[0].iov_len = sizeof(protocol_type);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);
#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\tTransmitting field characteristic ... ";
#endif

  uint64_t q = FieldType::p;

  iov[0].iov_base = &q;
  iov[0].iov_len = sizeof(q);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\tReceiving protocol parameters... ";
#endif

  iov[0].iov_base = &this->set_;
  iov[0].iov_len = sizeof(this->set_);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\t\tM: " << this->set_.M << std::endl;
  std::cout << "\t\ttau: " << this->set_.tau << std::endl;
  std::cout << "\t\tN: " << this->set_.N << std::endl;
  std::cout << "\t\tn: " << this->set_.n << std::endl;
  std::cout << "\t\tm: " << this->set_.m << std::endl;
  std::cout << "\tReceiving protocol public data... ";
#endif

  iovec *iov2 = new iovec[set_.n + 1];

  a_.resize(set_.n);
  for (auto i = 0; i < set_.n; ++i) {
    a_[i].resize(set_.m);
    iov2[i].iov_base = a_[i].data();
    iov2[i].iov_len = a_[i].size() * sizeof(a_[i][0]);
  }
  t_.resize(set_.n);
  iov2[set_.n].iov_base = t_.data();
  iov2[set_.n].iov_len = t_.size() * sizeof(t_[0]);
  nread = readv(this->sock_, iov2, set_.n + 1);
  assert (nread == (int)iov2[0].iov_len * set_.n + (int)iov2[set_.n].iov_len);

  delete[] iov2;

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "Negotiating protocol parameters... done" << std::endl;
#endif

  return 0;
}

template<class FieldType>
bool CacVerifierParty<FieldType>::runOnline() {
  CacVerifierLogic<FieldType> v(set_, a_, t_);

  iovec *iov = new iovec[3000];
  ssize_t nwritten, nread;

  // ** Round 1 output **
  block h_gamma;
  iov[0].iov_base = &h_gamma;
  iov[0].iov_len = sizeof(h_gamma);
  nread = readv(this->sock_, iov, 1);
  assert (nread == (int)iov[0].iov_len);

  // ** Round 2 **
  std::vector<uint8_t> E;
  v.r2(h_gamma, E); // Run round 2
  iov[0].iov_base = E.data();
  iov[0].iov_len = E.size() * sizeof(E[0]);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  // ** Round 3 output **
  std::vector<block> seed(set_.tau), omegaN(set_.M - set_.tau);
  block h_pi;
  iov[0].iov_base = seed.data();
  iov[0].iov_len = seed.size() * sizeof(seed[0]);
  iov[1].iov_base = omegaN.data();
  iov[1].iov_len = omegaN.size() * sizeof(omegaN[0]);
  iov[2].iov_base = &h_pi;
  iov[2].iov_len = sizeof(h_pi);
  nread = readv(this->sock_, iov, 3);
  assert (nread == (int)(iov[0].iov_len + iov[1].iov_len + iov[2].iov_len));

  // ** Round 4 **
  block seed_ell;
  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
  iov[0].iov_base = &seed_ell;
  iov[0].iov_len = sizeof(seed_ell);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  // ** Round 5 output **
  block h_psi;
  iov[0].iov_base = &h_psi;
  iov[0].iov_len = sizeof(h_psi);
  nread = readv(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  // ** Round 6 **
  std::vector<int> i_bar;
  v.r6(h_psi, i_bar); // Run round 6
  iov[0].iov_base = i_bar.data();
  iov[0].iov_len = i_bar.size() * sizeof(i_bar[0]);
  nwritten = writev(this->sock_, iov, 1);
  assert (nwritten == (int)iov[0].iov_len);

  // ** Round 7 output **
  block seed_e_bar;
  std::vector<std::vector<block>> seed_tree(set_.M - set_.tau);
  std::vector<block> gamma_i_bar(set_.M - set_.tau);
  std::vector<std::vector<FieldType>> alpha_i_bar(set_.M - set_.tau), b_square(set_.M - set_.tau), s(set_.M - set_.tau);
  std::vector<FieldType> o_i_bar(set_.M - set_.tau);
  auto iov_id = 0;
  iov[iov_id].iov_base = &seed_e_bar;
  iov[iov_id++].iov_len = sizeof(seed_e_bar);
  for (auto i = 0; i < set_.M - set_.tau; ++i) {
    seed_tree[i].resize(log(set_.N)/log(2));
    iov[iov_id].iov_base = seed_tree[i].data();
    iov[iov_id++].iov_len = seed_tree[i].size() * sizeof(seed_tree[i][0]);
  }
  iov[iov_id].iov_base = gamma_i_bar.data();
  iov[iov_id++].iov_len = gamma_i_bar.size() * sizeof(gamma_i_bar[0]);
  for (auto i = 0; i < set_.M - set_.tau; ++i) {
    alpha_i_bar[i].resize(set_.m);
    iov[iov_id].iov_base = alpha_i_bar[i].data();
    iov[iov_id++].iov_len = alpha_i_bar[i].size() * sizeof(alpha_i_bar[i][0]);
  }
  iov[iov_id].iov_base = o_i_bar.data();
  iov[iov_id++].iov_len = o_i_bar.size() * sizeof(o_i_bar[0]);
  int i_id = 0;
  for (auto e = 0, e_id = 0; e < set_.M; ++e) {
    if (E[e])
      continue;

    if (v.verifiers_[e]->i_bar_ != set_.N - 1) {
      b_square[e_id].resize(set_.m);
      iov[iov_id].iov_base = b_square[e_id].data();
      iov[iov_id++].iov_len = b_square[e_id].size() * sizeof(b_square[e_id][0]);
      s[e_id].resize(set_.m);
      iov[iov_id].iov_base = s[e_id].data();
      iov[iov_id++].iov_len = s[e_id].size() * sizeof(s[e_id][0]);

      i_id += 2;
    }

    e_id++;
  }
  nread = readv(this->sock_, iov, 1 + (set_.M - set_.tau) + 1 + (set_.M - set_.tau) + 1 + i_id);
  assert (nread == (int)(iov[0].iov_len + iov[1].iov_len * (set_.M - set_.tau)  + iov[set_.M - set_.tau + 1].iov_len +
                         iov[set_.M - set_.tau + 2].iov_len * (set_.M - set_.tau) + iov[2 * (set_.M - set_.tau) + 2].iov_len +
                         iov[2 * (set_.M - set_.tau) + 3].iov_len * i_id));

  bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8

  if (flag)
    std::cout << "Proof accepted" << std::endl;
  else
    std::cout << "Proof rejected" << std::endl;

  delete[] iov;

  return true;
}


}


#endif // __LZKP_CAC_VERIFIER_PARTY_H_FILE__