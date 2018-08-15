//
// Created by roee on 7/26/18.
//

#ifndef LZKP_VERIFIER_H
#define LZKP_VERIFIER_H

#include <vector>

#include "settings.h"

namespace lzkp {

class Verifier {
public:
  Verifier(const Settings &s);

  const std::vector<bool> &r2();

private:
  const int M;
  const int N;
  const int q;
  const int m;
  const int t;

  std::vector<bool> E_;
};
}

#endif //LZKP_VERIFIER_H
