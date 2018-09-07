#include <iostream>


#include "factory.h"
#include "party.h"

//#include "settings.h"
//#include "seedtree.h"
//#include "Mersenne.h"

//#define PORT 20000

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;
using namespace lzkp;

/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
void option_dependency(const po::variables_map& vm,
                       const char* for_what, const char* required_option)
{
  if (vm.count(for_what) && !vm[for_what].defaulted())
    if (vm.count(required_option) == 0 || vm[required_option].defaulted())
      throw logic_error(string("Option '") + for_what
                        + "' requires option '" + required_option + "'.");
}

/* Function used to check that of 'for_what' is not specified, then
   'required_option' is specified. */
void option_negative_dependency(const po::variables_map& vm,
                       const char* for_what, const char* required_option)
{
  if (!vm.count(for_what) || vm[for_what].defaulted())
    if (vm.count(required_option) == 0 || vm[required_option].defaulted())
      throw logic_error(string("Lack of option '") + for_what
                        + "' requires option '" + required_option + "'.");
}

int main(int argc, char *argv[]) {
  int protocol_type;
  bool is_prover = false;

  int q, M, tau, N, n, m;
  bool is_accepted = false;

  string ip;
  int port;
//  uint64_t seed;

  po::variables_map vm;

  try {
    po::options_description generic("Generic options");
    generic.add_options()
      ("help,h", "print usage message")
      ("protocol_type,p", po::value<int>(&protocol_type)->required(), "protocol type (0=Cut-and-Choose, 1=Sacrificing)")
      ("is_prover,i", po::bool_switch(&is_prover), "running the prover?")
      ;

    po::options_description network("Network options");
    network.add_options()
      ("ip", po::value<string>(&ip), "other party IP (required only for verifier)")
      ("port", po::value<int>(&port)->required(), "port")
      ;

    po::options_description parameter("Parameter options");
    parameter.add_options()
      ("q,q", po::value<int>(&q), "field type (15=15BitField, 31=MersenneInt, 61=MersenneLong)")
      ("M,M", po::value<int>(&M), "number of MPCs 'in the head'")
      ("tau,t", po::value<int>(&tau), "number of games to open")
      ("N,N", po::value<int>(&N), "number of parties in each MPC")
      ("n,n", po::value<int>(&n), "rows in matrix")
      ("m,m", po::value<int>(&m), "columns in matrix")
      ("is_accepted,a", po::bool_switch(&is_accepted), "should proof be accepted?")
      ;

    po::options_description cmdline_options;

    cmdline_options.add(generic).add(network).add(parameter);

    store(parse_command_line(argc, argv, cmdline_options), vm);

    if (vm.count("help")) {
      cout << cmdline_options << endl;

      return 0;
    }

    option_dependency(vm, "is_prover", "q");
    option_dependency(vm, "is_prover", "M");
    option_dependency(vm, "is_prover", "tau");
    option_dependency(vm, "is_prover", "N");
    option_dependency(vm, "is_prover", "n");
    option_dependency(vm, "is_prover", "m");

    option_negative_dependency(vm, "is_prover", "ip");

    po::notify(vm);

    if (protocol_type != 0 && protocol_type != 1)
      throw po::validation_error(po::validation_error::invalid_option_value, "protocol_type");

    if (vm.count("q") && !(q == 15 || q == 31 || q == 61))
      throw po::validation_error(po::validation_error::invalid_option_value, "q");

    if (vm.count("M") && M < 2)
      throw po::validation_error(po::validation_error::invalid_option_value, "M");

    if (vm.count("tau") && !(tau > 0 && tau <= M))
      throw po::validation_error(po::validation_error::invalid_option_value, "tau");

    if (vm.count("N") && !(N > 0))
      throw po::validation_error(po::validation_error::invalid_option_value, "N");

    if (vm.count("n") && !(n > 0))
      throw po::validation_error(po::validation_error::invalid_option_value, "n");

    if (vm.count("m") && !(m > n))
      throw po::validation_error(po::validation_error::invalid_option_value, "m");

  }
  catch(exception &e) {
    cout << e.what() << endl;

    return 1;
  }

#ifdef DEBUG
  std::cout << "Protocol type: " << protocol_type << std::endl;
  std::cout << "IsProver? " << is_prover << std::endl;
  std::cout << "IP: " << ip << std::endl;
  std::cout << "Port: " << port << std::endl;
  std::cout << "q: " << q << std::endl;
#endif

  Party *party = Factory::create(protocol_type, is_prover, q);

  if (party->init(argc, argv)) {
#ifdef DEBUG
    std::cout << "protocol initialization error" << std::endl;
#endif

    delete party;

    return 1;
  }

  // party->runOnline();

  delete party;

  return 0;

  if (argc != 10) {
    cout << "usage: LZKP party_id M tau N n m seed accepted_proof? party_ip" << endl;

    return 0;
  }

//  int is_prover = !atoi(argv[1]);

//  int M       = atoi(argv[2]);
//  int tau     = atoi(argv[3]);
//  int N       = atoi(argv[4]);
//  int n       = atoi(argv[5]);
//  int m       = atoi(argv[6]);
//  block seed;
//  seed.halves[0] = 0;
//  seed.halves[1] = atoi(argv[7]);
//  int acc    = atoi(argv[8]);  // Should the proof be accepted?
//  string ip   = argv[9];


//  if (is_prover)
//    cout << "Starting prover" << endl;
//  else
//    cout << "Starting verifier" << endl;
//
//  Settings set(M, tau, N, ZpMersenneLongElement::p, n, m); // need to swap n,m
//
//  std::vector<std::vector<ZpMersenneLongElement>> a;
//  std::vector<ZpMersenneLongElement> t, secret;
//
//  }

//  if (!acc) {
//    for (auto mm = 0; mm < m; ++mm) {
//      secret[mm] += ZpMersenneLongElement(prng.get<block>().bytes[0] % 2);
//    }
//  }

  assert(M > tau);
  assert(N > 0);
  assert(n > 0);
  assert(m > n);
  assert(tau > 0);

  cout << "Creating communication channel... (" << ip << ")" << endl;

//

//    CacProverLogic<ZpMersenneLongElement> p(set, a, t, secret);
//
//    iovec *iov = new iovec[3000];
//    ssize_t nwritten, nread;
//
//    // ** Round 1 **
//    block h_gamma;
//    p.r1(h_gamma); // Run round 1
////    cout << "h_gamma " << h_gamma.halves[0] << endl;
//    iov[0].iov_base = &h_gamma;
//    iov[0].iov_len = sizeof(h_gamma);
//    nwritten = writev(new_socket, iov, 1);
//    assert (nwritten == (int)iov[0].iov_len);
//
//    // ** Round 2 output **
//    std::vector<uint8_t> E(M);
//    iov[0].iov_base = E.data();
//    iov[0].iov_len = E.size() * sizeof(E[0]);
//    nread = readv(new_socket, iov, 1);
//    assert (nread == (int)iov[0].iov_len);
//
//    // ** Round 3 **
//    std::vector<block> seed, omegaN;
//    block h_pi;
//    p.r3(E, seed, omegaN, h_pi); // Run round 3
//    iov[0].iov_base = seed.data();
//    iov[0].iov_len = seed.size() * sizeof(seed[0]);
//    iov[1].iov_base = omegaN.data();
//    iov[1].iov_len = omegaN.size() * sizeof(omegaN[0]);
//    iov[2].iov_base = &h_pi;
//    iov[2].iov_len = sizeof(h_pi);
//    nwritten = writev(new_socket, iov, 3);
//    assert (nwritten == (int)(iov[0].iov_len + iov[1].iov_len + iov[2].iov_len));
//
//    // ** Round 4 output**
//    block seed_ell;
//    iov[0].iov_base = &seed_ell;
//    iov[0].iov_len = sizeof(seed_ell);
//    nread = readv(new_socket, iov, 1);
//    assert (nread == (int)iov[0].iov_len);
//
//    // ** Round 5 **
//    block h_psi;
//    p.r5(seed_ell, h_psi); // Run round 5
//    iov[0].iov_base = &h_psi;
//    iov[0].iov_len = sizeof(h_psi);
//    nwritten = writev(new_socket, iov, 1);
//    assert (nwritten == (int)iov[0].iov_len);
//
//    // ** Round 6 output **
//    std::vector<int> i_bar(M - tau);
//    iov[0].iov_base = i_bar.data();
//    iov[0].iov_len = i_bar.size() * sizeof(i_bar[0]);
//    nread = readv(new_socket, iov, 1);
//    assert (nread == (int)iov[0].iov_len);
//
//    // ** Round 7 **
//    block seed_e_bar;
//    std::vector<std::vector<block>> seed_tree;
//    std::vector<block> gamma_i_bar;
//    std::vector<std::vector<ZpMersenneLongElement>> alpha_i_bar, b_square, s;
//    std::vector<ZpMersenneLongElement> o_i_bar;
//    p.r7(i_bar, seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 7
//    auto iov_id = 0;
//    iov[iov_id].iov_base = &seed_e_bar;
//    iov[iov_id++].iov_len = sizeof(seed_e_bar);
//    for (auto i = 0; i < M - tau; ++i) {
//      iov[iov_id].iov_base = seed_tree[i].data();
//      iov[iov_id++].iov_len = seed_tree[i].size() * sizeof(seed_tree[i][0]);
//    }
//    iov[iov_id].iov_base = gamma_i_bar.data();
//    iov[iov_id++].iov_len = gamma_i_bar.size() * sizeof(gamma_i_bar[0]);
//    for (auto i = 0; i < M - tau; ++i) {
//      iov[iov_id].iov_base = alpha_i_bar[i].data();
//      iov[iov_id++].iov_len = alpha_i_bar[i].size() * sizeof(alpha_i_bar[i][0]);
//    }
//    iov[iov_id].iov_base = o_i_bar.data();
//    iov[iov_id++].iov_len = o_i_bar.size() * sizeof(o_i_bar[0]);
//    int i_id = 0;
//    for (auto e = 0, e_id = 0; e < M; ++e) {
//      if (E[e])
//        continue;
//
//      if (p.provers_[e]->i_bar_ != N - 1) {
//        iov[iov_id].iov_base = b_square[e_id].data();
//        iov[iov_id++].iov_len = b_square[e_id].size() * sizeof(b_square[e_id][0]);
//        iov[iov_id].iov_base = s[e_id].data();
//        iov[iov_id++].iov_len = s[e_id].size() * sizeof(s[e_id][0]);
//
//        i_id += 2;
//      }
//
//      e_id++;
//    }
//    nwritten = writev(new_socket, iov, 1 + (M - tau) + 1 + (M - tau) + 1 + i_id);
//    assert (nwritten == (int)(iov[0].iov_len + iov[1].iov_len * (M - tau) + iov[M - tau + 1].iov_len +
//                              iov[M - tau + 2].iov_len * (M - tau) + iov[2 * (M - tau) + 2].iov_len +
//                              iov[2 * (M - tau) + 3].iov_len * i_id));
//
//    delete[] iov;
//    close(new_socket);
//    cout << "Done!" << endl;
//  }
//  else {
//
//
//    CacVerifierLogic<ZpMersenneLongElement> v(set, a, t);
//
//    iovec *iov = new iovec[3000];
//    ssize_t nwritten, nread;
//
//    // ** Round 1 output **
//    block h_gamma;
//    iov[0].iov_base = &h_gamma;
//    iov[0].iov_len = sizeof(h_gamma);
//    nread = readv(sock, iov, 1);
//    assert (nread == (int)iov[0].iov_len);
////    cout << "h_gamma " << h_gamma.halves[0] << endl;
//
//    // ** Round 2 **
//    std::vector<uint8_t> E;
//    v.r2(h_gamma, E); // Run round 2
//    iov[0].iov_base = E.data();
//    iov[0].iov_len = E.size() * sizeof(E[0]);
//    nwritten = writev(sock, iov, 1);
//    assert (nwritten == (int)iov[0].iov_len);
//
//    // ** Round 3 output **
//    std::vector<block> seed(tau), omegaN(M - tau);
//    block h_pi;
//    iov[0].iov_base = seed.data();
//    iov[0].iov_len = seed.size() * sizeof(seed[0]);
//    iov[1].iov_base = omegaN.data();
//    iov[1].iov_len = omegaN.size() * sizeof(omegaN[0]);
//    iov[2].iov_base = &h_pi;
//    iov[2].iov_len = sizeof(h_pi);
//    nread = readv(sock, iov, 3);
//    assert (nread == (int)(iov[0].iov_len + iov[1].iov_len + iov[2].iov_len));
//
//    // ** Round 4 **
//    block seed_ell;
//    v.r4(seed, omegaN, h_pi, seed_ell); // Run round 4
//    iov[0].iov_base = &seed_ell;
//    iov[0].iov_len = sizeof(seed_ell);
//    nwritten = writev(sock, iov, 1);
//    assert (nwritten == (int)iov[0].iov_len);
//
//    // ** Round 5 output **
//    block h_psi;
//    iov[0].iov_base = &h_psi;
//    iov[0].iov_len = sizeof(h_psi);
//    nread = readv(sock, iov, 1);
//    assert (nwritten == (int)iov[0].iov_len);
//
//    // ** Round 6 **
//    std::vector<int> i_bar;
//    v.r6(h_psi, i_bar); // Run round 6
//    iov[0].iov_base = i_bar.data();
//    iov[0].iov_len = i_bar.size() * sizeof(i_bar[0]);
//    nwritten = writev(sock, iov, 1);
//    assert (nwritten == (int)iov[0].iov_len);
//
//    // ** Round 7 output **
//    block seed_e_bar;
//    std::vector<std::vector<block>> seed_tree(M - tau);
//    std::vector<block> gamma_i_bar(M - tau);
//    std::vector<std::vector<ZpMersenneLongElement>> alpha_i_bar(M - tau), b_square(M - tau), s(M - tau);
//    std::vector<ZpMersenneLongElement> o_i_bar(M - tau);
//    auto iov_id = 0;
//    iov[iov_id].iov_base = &seed_e_bar;
//    iov[iov_id++].iov_len = sizeof(seed_e_bar);
//    for (auto i = 0; i < M - tau; ++i) {
//      seed_tree[i].resize(log(N)/log(2));
//      iov[iov_id].iov_base = seed_tree[i].data();
//      iov[iov_id++].iov_len = seed_tree[i].size() * sizeof(seed_tree[i][0]);
//    }
//    iov[iov_id].iov_base = gamma_i_bar.data();
//    iov[iov_id++].iov_len = gamma_i_bar.size() * sizeof(gamma_i_bar[0]);
//    for (auto i = 0; i < M - tau; ++i) {
//      alpha_i_bar[i].resize(m);
//      iov[iov_id].iov_base = alpha_i_bar[i].data();
//      iov[iov_id++].iov_len = alpha_i_bar[i].size() * sizeof(alpha_i_bar[i][0]);
//    }
//    iov[iov_id].iov_base = o_i_bar.data();
//    iov[iov_id++].iov_len = o_i_bar.size() * sizeof(o_i_bar[0]);
//    int i_id = 0;
//    for (auto e = 0, e_id = 0; e < M; ++e) {
//      if (E[e])
//        continue;
//
//      if (v.verifiers_[e]->i_bar_ != N - 1) {
//        b_square[e_id].resize(m);
//        iov[iov_id].iov_base = b_square[e_id].data();
//        iov[iov_id++].iov_len = b_square[e_id].size() * sizeof(b_square[e_id][0]);
//        s[e_id].resize(m);
//        iov[iov_id].iov_base = s[e_id].data();
//        iov[iov_id++].iov_len = s[e_id].size() * sizeof(s[e_id][0]);
//
//        i_id += 2;
//      }
//
//      e_id++;
//    }
//    nread = readv(sock, iov, 1 + (M - tau) + 1 + (M - tau) + 1 + i_id);
//    assert (nread == (int)(iov[0].iov_len + iov[1].iov_len * (M - tau)  + iov[M - tau + 1].iov_len +
//                           iov[M - tau + 2].iov_len * (M - tau) + iov[2 * (M - tau) + 2].iov_len +
//                           iov[2 * (M - tau) + 3].iov_len * i_id));
//
//    bool flag = v.r8(seed_e_bar, seed_tree, gamma_i_bar, alpha_i_bar, o_i_bar, b_square, s); // Run round 8
//
//    if (flag)
//      cout << "Proof accepted" << endl;
//    else
//      cout << "Proof rejected" << endl;
//
//    delete[] iov;
//    close(sock);
//    cout << "Done!" << endl;
//  }

  return 0;
}