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
class SacProverParty : public ProverParty {
public:
  SacProverParty() : ProverParty() {
    debug("Constructing SacProverParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl);
  }
  ~SacProverParty() {
    debug("Destructing SacProverParty<" << boost::typeindex::type_id<FieldType>().pretty_name() << ">" << std::endl);

    free2Darray<FieldType>(a_);
    free1Darray<FieldType>(t_);
    free1Darray<FieldType>(secret_);
  }

  virtual int init(int argc, const char* const argv[]);
  virtual bool runOnline();

  static const int PROTOCOL_TYPE = 1;

protected:
  virtual int parseArguments(int argc, const char* const argv[]);
  virtual int negotiateParameters();
  virtual int generateData();

  // Public known values
  FieldType **a_;
  FieldType *t_;

  // Prover's secret
  FieldType *secret_;

  Parameters par_;
  int M;
  int N;
  int n;
  int m;
  bool is_accepted = false;

  int x_;
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

  this->sync();

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
      ("multi_threaded,x", po::value<int>(&x_)->default_value(1), "number of threads")
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

  if (x_ != 1)
    debug("\tMulti-threading enabled (# " << x_ << " threads)" << std::endl);
  else
    debug("\tMultit-hreading disabled" << std::endl);

  return 0;
}

template<class FieldType>
int SacProverParty<FieldType>::negotiateParameters() {
  debug("Negotiating protocol parameters..." << std::endl);

  iovec iov[2];

  int protocol_type;

  iov[0].iov_base = &protocol_type;
  iov[0].iov_len = sizeof(protocol_type);
  this->readvWrapper(iov, 1, iov[0].iov_len);

  debug("\tValidating protocol type... ");

  assert (protocol_type == SacProverParty::PROTOCOL_TYPE);

  debug("done" << std::endl);

  uint64_t q;

  iov[0].iov_base = &q;
  iov[0].iov_len = sizeof(q);
  this->readvWrapper(iov, 1, iov[0].iov_len);

  debug("\tValidating field... ");

  assert (q == (uint64_t)FieldType::p);

  par_ = Parameters(M, 0, N, n, m);

  debug("done" << std::endl);
  debug("\tTransmitting protocol parameters... ");

  iov[0].iov_base = &par_;
  iov[0].iov_len = sizeof(par_);
  this->writevWrapper(iov, 1, iov[0].iov_len);

  debug("done" << std::endl);

  this->generateData(); // Generate protocol data

  debug("\tTransmitting protocol public data... ");

  iov[0].iov_base = a_[0];
  iov[0].iov_len = n * m * sizeof(FieldType);
  iov[1].iov_base = t_;
  iov[1].iov_len = n * sizeof(FieldType);
  this->writevWrapper(iov, 2, iov[0].iov_len + iov[1].iov_len);

  debug("done" << std::endl);
  debug("Negotiating protocol parameters... done" << std::endl);

  return 0;
}

template<class FieldType>
int SacProverParty<FieldType>::generateData() {
  osuCrypto::PRNG prng(osuCrypto::sysRandomSeed());

  a_ = allocate2D<FieldType>(n, m);
  t_ = allocate1D<FieldType>(n);
  secret_ = allocate1D<FieldType>(m);

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
  SacProverLogic<FieldType> p(par_, a_, t_, secret_, x_);

  iovec iov[9];

  // ** Round 1 **
  block h_gamma;
  debug("\tExecuting round #1... ");
  startComputationClock();
  p.r1(h_gamma); // Run round 1
  stopComputationClock();
  debug("done" << std::endl);
  iov[0].iov_base = &h_gamma;
  iov[0].iov_len = sizeof(h_gamma);
  debug("\tSending output of round #1... ");
  this->writevWrapper(iov, 1, iov[0].iov_len);
  debug("done" << std::endl);

  // ** Round 2 output **
  block seed_ell;
  iov[0].iov_base = &seed_ell;
  iov[0].iov_len = sizeof(seed_ell);
  debug("\tReceiving output of round #2... ");
  this->readvWrapper(iov, 1, iov[0].iov_len);
  debug("done" << std::endl);

  // ** Round 3 **
  block h_pi, h_psi, h_theta;
  debug("\tExecuting round #3... ");
  startComputationClock();
  p.r3(seed_ell, h_pi, h_psi, h_theta); // Run round 3
  stopComputationClock();
  time_eq_1 = p.time_eq_1;
  tot_matrix_multiplication_time = p.tot_matrix_multiplication_time / 1000000;
  debug("done" << std::endl);
  iov[0].iov_base = &h_pi;
  iov[0].iov_len = sizeof(h_pi);
  iov[1].iov_base = &h_psi;
  iov[1].iov_len = sizeof(h_psi);
  iov[2].iov_base = &h_theta;
  iov[2].iov_len = sizeof(h_theta);
  debug("\tSending output of round #3... ");
  this->writevWrapper(iov, 3, iov[0].iov_len + iov[1].iov_len + iov[2].iov_len);
  debug("done" << std::endl);

  // ** Round 4 output **
  std::vector<int> i_bar(M);
  iov[0].iov_base = i_bar.data();
  iov[0].iov_len = i_bar.size() * sizeof(i_bar[0]);
  debug("\tReceiving output of round #4... ");
  this->readvWrapper(iov, 1, iov[0].iov_len);
  debug("done" << std::endl);

  // ** Round 5 **
  block seed_global;
  std::vector<std::vector<block>> seed_tree;
  std::vector<block> gamma_i_bar;
  std::vector<std::vector<FieldType>> alpha_i_bar, b_square, s, s_square;
  std::vector<FieldType> o_i_bar, v_i_bar;
  debug("\tExecuting round #5... ");
  startComputationClock();
  p.r5(i_bar, seed_global, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, v_i_bar, b_square, s, s_square); // Run round 5
  stopComputationClock();
  debug("done" << std::endl);

  int i_id = 0;
  for (auto e = 0; e < M; ++e) {
    if (p.provers_[e]->i_bar_ != N - 1) {
      i_id++;
    }
  }

  block **seed_tree_f = allocate2D<block>(M, (int)(log(N) / log(2)));
  FieldType **alpha_i_bar_f = allocate2D<FieldType>(M, m);
  FieldType **b_square_f = allocate2D<FieldType>(i_id, m);
  FieldType **s_f = allocate2D<FieldType>(i_id, m);
  FieldType **s_square_f = allocate2D<FieldType>(i_id, m);

  for (auto i = 0; i < M; ++i) {
    memcpy(seed_tree_f[i], seed_tree[i].data(), (int)(log(N) / log(2)) * sizeof(block));
    memcpy(alpha_i_bar_f[i], alpha_i_bar[i].data(), m * sizeof(FieldType));
  }

  i_id = 0;
  for (auto e = 0, e_id = 0; e < M; ++e) {
    if (p.provers_[e]->i_bar_ != N - 1) {
      memcpy(b_square_f[i_id], b_square[e].data(), m * sizeof(FieldType));
      memcpy(s_f[i_id], s[e].data(), m * sizeof(FieldType));
      memcpy(s_square_f[i_id], s_square[e].data(), m * sizeof(FieldType));

      i_id++;
    }

    e_id++;
  }

  iov[0].iov_base = &seed_global;
  iov[0].iov_len = sizeof(seed_global);
  iov[1].iov_base = seed_tree_f[0];
  iov[1].iov_len = M * (int)(log(N) / log(2)) * sizeof(block);
  iov[2].iov_base = gamma_i_bar.data();
  iov[2].iov_len = gamma_i_bar.size() * sizeof(gamma_i_bar[0]);
  iov[3].iov_base = alpha_i_bar_f[0];
  iov[3].iov_len = M * m * sizeof(FieldType);
  iov[4].iov_base = o_i_bar.data();
  iov[4].iov_len = o_i_bar.size() * sizeof(o_i_bar[0]);
  iov[5].iov_base = v_i_bar.data();
  iov[5].iov_len = v_i_bar.size() * sizeof(v_i_bar[0]);
  iov[6].iov_base = b_square_f[0];
  iov[6].iov_len = i_id * m * sizeof(FieldType);
  iov[7].iov_base = s_f[0];
  iov[7].iov_len = i_id * m * sizeof(FieldType);
  iov[8].iov_base = s_square_f[0];
  iov[8].iov_len = i_id * m * sizeof(FieldType);
  debug("\tSending output of round #5... ");
  this->writevWrapper(iov, 9, iov[0].iov_len + iov[1].iov_len + iov[2].iov_len + iov[3].iov_len +
                              iov[4].iov_len + iov[5].iov_len + iov[6].iov_len + iov[7].iov_len +
                              iov[8].iov_len);
  debug("done" << std::endl);

  bool flag;
  iov[0].iov_base = &flag;
  iov[0].iov_len = sizeof(flag);
  debug("\tReceiving protocol output... ");
  this->readvWrapper(iov, 1, iov[0].iov_len);
  debug("done" << std::endl);

  assert(flag == is_accepted);

  debug("Online phase... done" << std::endl);

  free2Darray<block>(seed_tree_f);
  free2Darray<FieldType>(alpha_i_bar_f);
  free2Darray<FieldType>(b_square_f);
  free2Darray<FieldType>(s_f);
  free2Darray<FieldType>(s_square_f);

  return true;
}


}


#endif // __LZKP_SAC_PROVER_PARTY_H_FILE__
