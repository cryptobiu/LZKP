#include <iostream>
#include <cstdlib>

#include "prover_wrapper.h"
#include "verifier_wrapper.h"
#include "settings.h"

using namespace std;
using namespace lzkp;

int main(int argc, char *argv[]) {
    if (argc != 10) {
        cout << "usage: LZKP party_id M N q n m tau accepted_proof?" << endl;
    }

    bool is_prover = !atoi(argv[1]);

    int M       = atoi(argv[2]);
    int N       = atoi(argv[3]);
    uint64_t q  = atoi(argv[4]);
    int n       = atoi(argv[5]);
    int m       = atoi(argv[6]);
    int tau     = atoi(argv[7]);
    uint64_t s  = atoi(argv[8]);
    bool acc    = atoi(argv[9]);


    if (is_prover)
        cout << "Starting prover" << endl;
    else
        cout << "Starting verifier" << endl;

    Settings set(M, N, q, m, n, tau); // need to swap n,m

    NTL::ZZ_p::init(NTL::ZZ(q)); // *** CHECK HOW TO INIT 128 BIT ***

    NTL::Mat<NTL::ZZ_p> a;
    NTL::Vec<NTL::ZZ_p> t, secret;
    osuCrypto::PRNG prng(s);

    a.SetDims(n, m); // Fill with values
    t.SetLength(n); // Fill with values
    secret.SetLength(m);

    // Random matrix A
    for (auto nn = 0; nn < n; ++nn) {
        for (auto mm = 0; mm < m; ++mm) {
            a[nn][mm] = NTL::ZZ_p(prng.get<block>().halves[0]);
        }
    }

    // Random vector secret
    for (auto mm = 0; mm < m; ++mm) {
        secret[mm] = prng.get<block>().bytes[0] % 2;
    }

    // Calculate t
    for (auto nn = 0; nn < n; ++nn) {
        t[nn] = 0;

        for (auto mm = 0; mm < m; ++mm) {
            t[nn] += a[nn][mm] * secret[mm];
        }
    }

    return 0;
}