#ifndef _bueht_rotate
#define _bueht_rotate

#include "bueht_parameters.h"
#include "basis_set.h"
#include "atom.h"

/*
  James McNeely
*/

namespace BUEHT
{
  void Overlap ( const BUEHT::Atom & atom1, const BUEHT::Atom & atom2,
                    const BUEHT::BasisSet & basis1,
                    const BUEHT::BasisSet & basis2,
                    double* overlap_matirx )
  {
    int num_cart1, num_cart2;
    for( unsigned int i = 0; i < basis1.GetNumShells(); i++ )
    {
      for ( unsigned int j = 0; j < basis2.GetNumShells(); j++ )
      {
        num_cart1 = BUEHT::CartesianExpansionLength(basis1.GetL(i),0);
        num_cart2 = BUEHT::CartesianExpansionLength(basis2.GetL(j),0);
        std::cout << "L1: " << num_cart1 << "\n";
        std::cout << "L2: " << num_cart2 << "\n";
      }
    }
  }
} // End of Namespace BUEHT

#endif
