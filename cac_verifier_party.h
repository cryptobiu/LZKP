//
// Created by lzkp on 9/7/18.
//

#ifndef __LZKP_CAC_VERIFIER_PARTY_H_FILE__
#define __LZKP_CAC_VERIFIER_PARTY_H_FILE__


#include "verifier_party.h"

#include "parameters.h"
#include "cac_verifier_logic.h"

#include <cstring>
namespace lzkp {


template <class FieldType>
class CacVerifierParty : public VerifierParty {
public:
  CacVerifierParty() : VerifierParty(), a_(nullptr), t_(nullptr) {
    debug("Constructing CacVerifierParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl);
  }

  ~CacVerifierParty() {
    debug("Destructing CacVerifierParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl);

    free2Darray<FieldType>(a_);
    free1Darray<FieldType>(t_);
  }

  virtual int init(int argc, const char* const argv[]);
  virtual bool runOnline();

  static const int PROTOCOL_TYPE = 0;

protected:
  virtual int parseArguments(int argc, const char* const argv[]);
  virtual int negotiateParameters();

  // Public known values
  FieldType **a_;
  FieldType *t_;

  Parameters par_;

  int x_;
};

template<class FieldType>
int CacVerifierParty<FieldType>::init(int argc, const char* const argv[]) {
  parseArguments(argc, argv);
  int com = this->initCommunication();

  if (com != 0)
    std::cerr << "Cannot connect to prover (error " << com << ")" << std::endl;

  int par = this->negotiateParameters();

  if (par != 0)
    std::cerr << "Parameters negotiation failed (error " << par << ")" << std::endl;

  this->sync();

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

  po::options_description performence("Performence options");
  performence.add_options()
      ("multi_threaded,x", po::value<int>(&x_)->default_value(1), "number of threads")
    ;

  po::options_description cmdline_options;

  cmdline_options.add(network).add(performence);

  po::store(po::command_line_parser(argc, argv).options(cmdline_options).allow_unregistered().run(), vm);

  notify(vm);

  debug("Parsed parameters: " << std::endl);
  debug("\tIP: " << this->ip_ << std::endl);
  debug("\tPort: " << this->port_ << std::endl);

  if (x_ != 1)
    debug("\tMulti-threading enabled (# " << x_ << " threads)" << std::endl);
  else
    debug("\tMultit-hreading disabled" << std::endl);

  return 0;
}

template<class FieldType>
int CacVerifierParty<FieldType>::negotiateParameters() {
  debug("Negotiating protocol parameters..." << std::endl);

  iovec iov[2];

  int protocol_type = CacVerifierParty::PROTOCOL_TYPE;

  debug("\tTransmitting protocol type... ");

  iov[0].iov_base = &protocol_type;
  iov[0].iov_len = sizeof(protocol_type);
  this->writevWrapper(iov, 1, iov[0].iov_len);

  debug("done" << std::endl);
  debug("\tTransmitting field characteristic ... ");

  uint64_t q = FieldType::p;

  iov[0].iov_base = &q;
  iov[0].iov_len = sizeof(q);
  this->writevWrapper(iov, 1, iov[0].iov_len);

  debug("done" << std::endl);
  debug("\tReceiving protocol parameters... ");

  iov[0].iov_base = &this->par_;
  iov[0].iov_len = sizeof(this->par_);
  this->readvWrapper(iov, 1, iov[0].iov_len);

  debug("done" << std::endl);
  debug("\t\tM: " << this->par_.M << std::endl);
  debug("\t\ttau: " << this->par_.tau << std::endl);
  debug("\t\tN: " << this->par_.N << std::endl);
  debug("\t\tn: " << this->par_.n << std::endl);
  debug("\t\tm: " << this->par_.m << std::endl);
  debug("\tReceiving protocol public data... ");

  a_ = allocate2D<FieldType>(par_.n, par_.m);
  t_ = allocate1D<FieldType>(par_.n);

  iov[0].iov_base = a_[0];
  iov[0].iov_len = par_.n * par_.m * sizeof(FieldType);
  iov[1].iov_base = t_;
  iov[1].iov_len = par_.n * sizeof(FieldType);
  this->readvWrapper(iov, 2, iov[0].iov_len + iov[1].iov_len);

  debug("done" << std::endl);
  debug("Negotiating protocol parameters... done" << std::endl);

  return 0;
}

template<class FieldType>
bool CacVerifierParty<FieldType>::runOnline() {
  CacVerifierLogic<FieldType> v(par_, a_, t_, x_);

  iovec iov[7];

  debug("Online phase:" << std::endl);

  // ** Round 1 output **
  block h_gamma;
  iov[0].iov_base = &h_gamma;
  iov[0].iov_len = sizeof(h_gamma);
  debug("\tReceiving output of round #1... ");
  this->readvWrapper(iov, 1, iov[0].iov_len);
  debug("done" << std::endl);

  // ** Round 2 **
  std::vector<uint8_t> E;
  debug("\tExecuting round #2... ");
  startComputationClock();
  v.r2(h_gamma, E); // Run round 2
  stopComputationClock();
  debug("done" << std::endl);
  iov[0].iov_base = E.data();
  iov[0].iov_len = E.size() * sizeof(E[0]);
  debug("\tSending output of round #2... ");
  this->writevWrapper(iov, 1, iov[0].iov_len);
  debug("done" << std::endl);

  // ** Round 3 output **
  std::vector<block> seed(par_.tau), omegaN(par_.M - par_.tau);
  block h_pi;
  iov[0].iov_base = seed.data();
  iov[0].iov_len = seed.size() * sizeof(seed[0]);
  iov[1].iov_base = omegaN.data();
  iov[1].iov_len = omegaN.size() * sizeof(omegaN[0]);
  iov[2].iov_base = &h_pi;
  iov[2].iov_len = sizeof(h_pi);
  debug("\tReceiving output of round #3... ");
  this->readvWrapper(iov, 3, iov[0].iov_len + iov[1].iov_len + iov[2].iov_len);
  debug("done" << std::endl);

//   ** Round 4 **
  block seed_ell;
  debug("\tExecuting round #4... ");
  startComputationClock();
  auto cut_and_choose_clock = std::chrono::high_resolution_clock::now();
  v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
  time_cut_and_choose_ = std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::high_resolution_clock::now() - cut_and_choose_clock).count();
  stopComputationClock();
  debug("done" << std::endl);
  iov[0].iov_base = &seed_ell;
  iov[0].iov_len = sizeof(seed_ell);
  debug("\tSending output of round #4... ");
  this->writevWrapper(iov, 1, iov[0].iov_len);
  debug("done" << std::endl);

  // ** Round 5 output **
  block h_psi;
  iov[0].iov_base = &h_psi;
  iov[0].iov_len = sizeof(h_psi);
  debug("\tReceiving output of round #5... ");
  this->readvWrapper(iov, 1, iov[0].iov_len);
  debug("done" << std::endl);

  // ** Round 6 **
  std::vector<int> i_bar;
  debug("\tExecuting round #6... ");
  startComputationClock();
  v.r6(h_psi, i_bar); // Run round 6
  stopComputationClock();
  debug("done" << std::endl);
  iov[0].iov_base = i_bar.data();
  iov[0].iov_len = i_bar.size() * sizeof(i_bar[0]);
  debug("\tSending output of round #6... ");
  this->writevWrapper(iov, 1, iov[0].iov_len);
  debug("done" << std::endl);

  // ** Round 7 output **
  block seed_e_bar;
  std::vector<std::vector<block>> seed_tree(par_.M - par_.tau);
  std::vector<block> gamma_i_bar(par_.M - par_.tau);
  std::vector<std::vector<FieldType>> alpha_i_bar(par_.M - par_.tau), b_square(par_.M - par_.tau), s(par_.M - par_.tau);
  std::vector<FieldType> o_i_bar(par_.M - par_.tau);

  int i_id = 0;
  for (auto e = 0; e < par_.M; ++e) {
    if (E[e])
      continue;

    if (v.verifiers_[e]->i_bar_ != par_.N - 1) {
      i_id++;
    }
  }

  block **seed_tree_f = allocate2D<block>(par_.M - par_.tau, (int)(log(par_.N) / log(2)));
  FieldType **alpha_i_bar_f = allocate2D<FieldType>(par_.M - par_.tau, par_.m);
  FieldType **b_square_f = allocate2D<FieldType>(i_id, par_.m);
  FieldType **s_f = allocate2D<FieldType>(i_id, par_.m);

  iov[0].iov_base = &seed_e_bar;
  iov[0].iov_len = sizeof(seed_e_bar);
  iov[1].iov_base = seed_tree_f[0];
  iov[1].iov_len = (par_.M - par_.tau) * (int)(log(par_.N) / log(2)) * sizeof(block);
  iov[2].iov_base = gamma_i_bar.data();
  iov[2].iov_len = gamma_i_bar.size() * sizeof(gamma_i_bar[0]);
  iov[3].iov_base = alpha_i_bar_f[0];
  iov[3].iov_len = (par_.M - par_.tau) * par_.m * sizeof(FieldType);
  iov[4].iov_base = o_i_bar.data();
  iov[4].iov_len = o_i_bar.size() * sizeof(o_i_bar[0]);
  iov[5].iov_base = b_square_f[0];
  iov[5].iov_len = i_id * par_.m * sizeof(FieldType);
  iov[6].iov_base = s_f[0];
  iov[6].iov_len = i_id * par_.m * sizeof(FieldType);
  debug("\tReceiving output of round #7... ");
  this->readvWrapper(iov, 7, iov[0].iov_len + iov[1].iov_len + iov[2].iov_len + iov[3].iov_len +
                             iov[4].iov_len + iov[5].iov_len + iov[6].iov_len);
  debug("done" << std::endl);

  for (auto i = 0; i < par_.M - par_.tau; ++i) {
    seed_tree[i].resize((int)(log(par_.N) / log(2)));
    alpha_i_bar[i].resize(par_.m);

    memcpy(seed_tree[i].data(), seed_tree_f[i], (int)(log(par_.N) / log(2)) * sizeof(block));
    memcpy(alpha_i_bar[i].data(), alpha_i_bar_f[i], par_.m * sizeof(FieldType));
  }

  i_id = 0;
  for (auto e = 0, e_id = 0; e < par_.M; ++e) {
    if (E[e])
      continue;

    if (v.verifiers_[e]->i_bar_ != par_.N - 1) {
      b_square[e_id].resize(par_.m);
      s[e_id].resize(par_.m);

      memcpy(b_square[e_id].data(), b_square_f[i_id], par_.m * sizeof(FieldType));
      memcpy(s[e_id].data(), s_f[i_id], par_.m * sizeof(FieldType));

      i_id++;
    }

    e_id++;
  }

  startComputationClock();
  this->is_accepted_ = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
  time_eq_1 = v.time_eq_1;
  tot_matrix_multiplication_time = v.tot_matrix_multiplication_time;
  stopComputationClock();

  debug(std::endl);
  if (this->is_accepted_)
          debug("\tProof accepted" << std::endl);
  else
          debug("\tProof rejected" << std::endl);

  iov[0].iov_base = &this->is_accepted_;
  iov[0].iov_len = sizeof(this->is_accepted_);
  debug("\tSending protocol output... ");
  this->writevWrapper(iov, 1, iov[0].iov_len);
  debug("done" << std::endl);

  debug("Online phase... done" << std::endl);

  free2Darray<block>(seed_tree_f);
  free2Darray<FieldType>(alpha_i_bar_f);
  free2Darray<FieldType>(b_square_f);
  free2Darray<FieldType>(s_f);

  return true;
}


}


#endif // __LZKP_CAC_VERIFIER_PARTY_H_FILE__