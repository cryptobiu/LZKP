#include <iostream>
#include <chrono>

#include "factory.h"
#include "party.h"

#ifdef DEBUG
  #define debug(x) std::cerr << x
#else
  #define debug(x)
#endif

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

  int q, M, tau, N, n, m, x;
  bool is_accepted = false;

  string ip;
  int port;

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

    po::options_description performence("Performence options");
    performence.add_options()
      ("multi_threaded,x", po::value<int>(&x)->default_value(1), ("number of threads (limited to " + to_string(std::thread::hardware_concurrency()) + ")").c_str())
      ;

    po::options_description cmdline_options;

    cmdline_options.add(generic).add(network).add(performence).add(parameter);

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

    if (vm.count("q") && !(q == 15 || q == 31 || q == 59 || q == 61))
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

    if (vm.count("multi_threaded") && !((x > 0) && (x <= (int)std::thread::hardware_concurrency())))
      throw po::validation_error(po::validation_error::invalid_option_value, "x", to_string(x));

  }
  catch(exception &e) {
    cout << e.what() << endl;

    return 1;
  }

  debug("Protocol type: " << protocol_type << std::endl);
  debug("IsProver? " << is_prover << std::endl);
  debug("IP: " << ip << std::endl);
  debug("Port: " << port << std::endl);
  debug("q: " << q << std::endl);

  Party *party = Factory::create(protocol_type, is_prover, q);

  if (party == nullptr) {
    return 1;
  }

  if (party->init(argc, argv)) {
    debug("protocol initialization error" << std::endl);

    delete party;

    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();
  party->runOnline();
  auto stop = std::chrono::high_resolution_clock::now();
  auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  if (is_prover) {
    if (protocol_type == 0) // Cut-and-Choose
      std::cout << dur.count() - (party->RTT_ / 1000) << "," << party->tot_computation_time_ << ","
                << party->time_eq_1 << "," << party->tot_matrix_multiplication_time << ","
                << party->time_cut_and_choose_ << "," << party->RTT_ << std::endl;
    else
      std::cout << dur.count() - (party->RTT_ / 1000) << "," << party->tot_computation_time_ << ","
                << party->time_eq_1 << "," << party->tot_matrix_multiplication_time << "," << party->RTT_ << std::endl;
  }
  else {
    std::cout << dur.count() << "," << party->tot_computation_time_ << ","
              << party->time_eq_1 << "," << party->tot_matrix_multiplication_time << ","
              << dynamic_cast<VerifierParty*>(party)->isAccepted() << std::endl;
  }

  delete party;

  return 0;
}
