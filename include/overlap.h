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

  double Overlap_PRIM ( const double & exp1, const double & exp2, int i1, int j1, int k1,
                        int i2, int j2, int k2 )
  {
    return 1.0;
  }

  double Overlap_BF (const BUEHT::Atom & atom1, const BUEHT::Atom & atom2,
                     const double & dx, const double & dy, const double & dz,
                     const std::vector<double> & coefs1, const std::vector<double> & coefs2,
                     const std::vector<double> & exps1, const std::vector<double> & exps2,
                     int i1, int j1, int k1, int i2, int j2, int k2)
  {
    double sum = 0.0;
    for ( unsigned int m = 0; m < coefs1.size(); m++ )
    {
      for ( unsigned int n = 0; n < coefs2.size(); n++ )
      {
        double mu = exps1[m] * exps2[n] / ( exps1[m] + exps2[n] );
        double p = exps1[m] + exps2[n];
        double xp = ( exps1[m] * atom1.GetX() + exps2[n] * atom2.GetX() ) / p;
        double yp = ( exps1[m] * atom1.GetY() + exps2[n] * atom2.GetY() ) / p;
        double zp = ( exps1[m] * atom1.GetZ() + exps2[n] * atom2.GetZ() ) / p;
        sum += coefs1[m]*coefs2[n];//*
             //  Overlap_PRIM(mu,p,dx,dy,dz,i1,j1,k1,i2,j2,k2);
      }
    }
    return 1.0;
  }

  void Overlap ( const BUEHT::Atom & atom1, const BUEHT::Atom & atom2,
                    const BUEHT::BasisSet & basis1,
                    const BUEHT::BasisSet & basis2,
                    double* overlap_matirx )
  {
    int index = 0;
    int num_cart1, num_cart2;
    double sum = 0.0;
    double dx(atom2.GetX()-atom1.GetX());
    double dy(atom2.GetY()-atom1.GetY());
    double dz(atom2.GetZ()-atom1.GetZ());
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
                        BUEHT::Overlap_BF (atom1, atom2, dx, dy, dz,
                                           basis1.GetCoefficients(i), basis2.GetCoefficients(j),
                                           basis1.GetExponents(i),basis2.GetExponents(j),
                                           BUEHT::RealSphericalHarmonics
                                           [
                                             BUEHT::RealSphericalHarmonicsPtr
                                             [
                                               basis1.GetL(i) * ( basis1.GetL(i) + 1 ) - k
                                             ] + m
                                           ].x,
                                           BUEHT::RealSphericalHarmonics
                                           [
                                             BUEHT::RealSphericalHarmonicsPtr
                                             [
                                               basis1.GetL(i) * ( basis1.GetL(i) + 1 ) - k
                                             ] + m
                                           ].y,
                                           BUEHT::RealSphericalHarmonics
                                           [
                                             BUEHT::RealSphericalHarmonicsPtr
                                             [
                                               basis1.GetL(i) * ( basis1.GetL(i) + 1 ) - k
                                             ] + m
                                           ].z,
                                           BUEHT::RealSphericalHarmonics
                                           [
                                             BUEHT::RealSphericalHarmonicsPtr
                                             [
                                               basis2.GetL(j) * ( basis2.GetL(j) + 1 ) - l
                                             ] + n
                                           ].x,
                                           BUEHT::RealSphericalHarmonics
                                           [
                                             BUEHT::RealSphericalHarmonicsPtr
                                             [
                                               basis2.GetL(j) * ( basis2.GetL(j) + 1 ) - l
                                             ] + n
                                           ].y,
                                           BUEHT::RealSphericalHarmonics
                                           [
                                             BUEHT::RealSphericalHarmonicsPtr
                                             [
                                               basis2.GetL(j) * ( basis2.GetL(j) + 1 ) - l
                                             ] + n
                                           ].z)
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
