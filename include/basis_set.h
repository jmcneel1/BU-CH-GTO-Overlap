#ifndef _bueht_basis_set
#define _bueht_basis_set

#include "basis_function.h"
#include "atomic_properties.h"

/*
  James McNeely
*/

namespace BUEHT
{

class BasisSet : public BasisFunction
{

  public:
    BasisSet( const std::string & fname, int atomicnum )
    {
      std::string element = BUEHT::bueht_atomic_names[atomicnum-1];
      std::ifstream inFile(fname);
      if ( !inFile.good() )
      {
        std::cout << "Couldn't open basis set file! Exiting...\n";
        std::exit(1);
      }
      inFile.close();
    }
  private:
    std::vector<BUEHT::BasisFunction> myBasisSet;
};

}

#endif
