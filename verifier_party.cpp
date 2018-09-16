//
// Created by roee on 9/6/18.
//

#include "verifier_party.h"


using namespace lzkp;



int VerifierParty::initCommunication() {
    struct sockaddr_in serv_addr;
    int r;

    debug("Initializing communication channel..." << std::endl);
    debug("\tCreating socket... ");

    if ((sock_ = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
        debug("error (" << sock_ << ")" << std::endl);

        return -1;
    }

    std::memset(&serv_addr, '0', sizeof(serv_addr));

    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = inet_addr(ip_.c_str());
    serv_addr.sin_port = htons(port_);

    // Convert IPv4 and IPv6 addresses from text to binary form
    debug("done (" << sock_ << ")" << std::endl);
    debug("\tConverting IP (" << this->ip_.c_str() << ") from text to binary form... ");

    if ((r = inet_pton(AF_INET, this->ip_.c_str(), &serv_addr.sin_addr)) <= 0) {
        debug("error (" << r << ")" << std::endl);

        return -2;
    }

    debug("done" << std::endl);
    debug("\tConnecting to " << ip_ << ":" << port_ << "... ");

    if ((r = connect(sock_, (struct sockaddr *)&serv_addr, sizeof(serv_addr))) < 0) {
        debug("error (" << r << ")" << std::endl);

        return -3;
    }
    debug("done" << std::endl);
    debug("\tmeasuring RTT... ");
    auto n = this->sync();
    assert(n == 2);
    debug("done" << std::endl);


    debug("Initializing communication channel... done" << std::endl);

    return 0;
}

int VerifierParty::sync() {
  char dummy;
  int n = 0;

  n += read(sock_, &dummy, 1);
  n += write(sock_, &dummy, 1);

  return n;
}