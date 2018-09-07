//
// Created by lzkp on 9/7/18.
//

#ifndef __LZKP_CAC_VERIFIER_PARTY_H_FILE__
#define __LZKP_CAC_VERIFIER_PARTY_H_FILE__


#include "verifier_party.h"

#include "settings.h"


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
  virtual bool runOnline() { return true; }

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


}


#endif // __LZKP_CAC_VERIFIER_PARTY_H_FILE__
