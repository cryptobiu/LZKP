//
// Created by lzkp on 9/7/18.
//

#ifndef __LZKP_CAC_PROVER_PARTY_H_FILE__
#define __LZKP_CAC_PROVER_PARTY_H_FILE__


#include "prover_party.h"

#include <iostream>

#include "settings.h"


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
  virtual bool runOnline() { return true; }

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

  Settings set_;
  int M;
  int tau;
  int N;
  int n;
  int m;
  bool is_accepted = false;
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

  cmdline_options.add(network).add(parameter);

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

  set_ = Settings(M, tau, N, n, m);

#ifdef DEBUG
  std::cout << "done" << std::endl;
  std::cout << "\tTransmitting protocol parameters... ";
#endif

  iov[0].iov_base = &set_;
  iov[0].iov_len = sizeof(set_);
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

}


#endif // __LZKP_CAC_PROVER_PARTY_H_FILE__
