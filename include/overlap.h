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
    int index = 0;
    int num_cart1, num_cart2;
    double sum = 0.0;
    for( unsigned int i = 0; i < basis1.GetNumShells(); i++ )
    {
      for ( unsigned int j = 0; j < basis2.GetNumShells(); j++ )
      {
        for ( int k = -1*basis1.GetL(i); k <= basis1.GetL(i); k++ )
        {
          for ( int l = -1*basis2.GetL(j); l <= basis2.GetL(j); l++ )
          {
            sum = 0.0;
            for ( unsigned int m = 0; 
                  m < BUEHT::CartesianExpansionLength(basis1.GetL(i),k); m++ )
            {
              for ( unsigned int n = 0; 
                    n < BUEHT::CartesianExpansionLength(basis2.GetL(j),l); n++ )
              {
                sum += (BUEHT::RealSphericalHarmonics
                        [
                          BUEHT::RealSphericalHarmonicsPtr
                          [
                            basis1.GetL(i) * ( basis1.GetL(i) + 1 ) - k
                          ] + m 
                        ].factor *
                        BUEHT::RealSphericalHarmonics
                        [
                          BUEHT::RealSphericalHarmonicsPtr
                          [
                            basis2.GetL(j) * ( basis2.GetL(j) + 1 ) - l
                          ] + n 
                        ].factor *
                       );
              }
            }
            std::cout << sum << "\n";
          }
        }
      }
    }
  }
} // End of Namespace BUEHT

#endif
