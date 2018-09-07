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

  party->runOnline();

  delete party;

  return 0;




//  }
//  else {
//
//
//
//
//
//  }

  return 0;
}