#ifndef _bueht_basis_set
#define _bueht_basis_set

#include "basis_function.h"

/*
  James McNeely
*/

namespace BUEHT
{

class BasisSet : public BasisFunction
{

  public:
  private:
    std::vector<BUEHT::BasisFunction> myBasisSet;
};

}

#endif
