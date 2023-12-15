#ifndef _bueht_rotate
#define _bueht_rotate

#include "bueht_parameters.h"
#include "basis_set.h"
#include "atom.h"

/*

  James McNeely

  This header contains the necessary functions to calculate 
  the overlap integrals between a pair of atoms. All BFs centered
  on the give atoms are calculated.

  We have that X = sum(i=1..NP, ci*xi) where xi are spherical GTOs

  Each spherical GTO is then expanded in its Cartesian coordinates with the 
  bueht_parameters object RealSphericalHarmonics.

  The primitive normalization constants are stored in the BF itself.

  We thus have that 

  S12 = int[sum(i=1..NP){ci*xi}*sum(j=1..NP){cj*xj}]
      = int[sum(i=1..NP){ci*sum(k=1..NC1){ck*x^m1y^n1z^o1*Exp[-zeta1*r^2]}}*
            sum(j=1..NP){cj*sum(l=1..NC2){cl*x^m2y^n2z^o2*Exp[-zeta2*r^2]}}]

*/

namespace BUEHT
{

  /*
  */

  double Overlap_Recursion ( const double & s00, const double & p, 
                             const int & i, const int & j, 
                             const double & xpa, const double & xpb )
  {
    if ( i < 0 || j < 0 )
    {
      return 0.;
    }
    else if ( i == 0 && j == 0 )
    {
      return s00;
    }
    else
    {
      if ( i == 0 )
      {
        return xpb*Overlap_Recursion(s00,p,i,j-1,xpa,xpb)+(1.e0/(2.e0*p)*((double)j-1.e0)*
                                     Overlap_Recursion(s00,p,i,j-2,xpa,xpb));
      }
      else if ( j == 0 )
      {
        return xpa*Overlap_Recursion(s00,p,i-1,j,xpa,xpb)+(1.e0/(2.e0*p)*((double)j-1.e0)*
                                     Overlap_Recursion(s00,p,i,j-1,xpa,xpb));
      }
      else
      {
        return xpa*Overlap_Recursion(s00,p,i-1,j,xpa,xpb) + 
               1.e0/(2.e0*p)*(((double)i-1.e0)*Overlap_Recursion(s00,p,i-2,j,xpa,xpb) +
                              (double)j * Overlap_Recursion(s00,p,i-1,j-1,xpa,xpb));
      }
    }
  }

  double Overlap_PRIM ( const double & mu, const double & p, 
                        const double & dx, const double & dy, const double & dz,
                        const int & i1, const int & j1, const int & k1, 
                        const int & i2, const int & j2, const int & k2,
                        const double & xpa, const double & ypa, const double & zpa,
                        const double & xpb, const double & ypb, const double & zpb )
  {
    double xab_squared = dx*dx;
    double yab_squared = dy*dy;
    double zab_squared = dz*dz;
    double s00x = std::sqrt(bueht_pi/p)*std::exp(-mu*xab_squared);
    double s00y = std::sqrt(bueht_pi/p)*std::exp(-mu*yab_squared);
    double s00z = std::sqrt(bueht_pi/p)*std::exp(-mu*zab_squared);
    return Overlap_Recursion(s00x,p,i1,i2,xpa,xpb)*Overlap_Recursion(s00y,p,j1,j2,ypa,ypb)*
           Overlap_Recursion(s00z,p,k1,k2,zpa,zpb);
  }

  double Overlap_BF (const BUEHT::Atom & atom1, const BUEHT::Atom & atom2,
                     const double & dx, const double & dy, const double & dz,
                     const BUEHT::BasisFunction & bf1,
                     const BUEHT::BasisFunction & bf2,
                     int i1, int j1, int k1, int i2, int j2, int k2)
  {
    double sum = 0.0;
    double mu, p, xp, yp, zp, xpa, ypa, zpa, xpb, ypb, zpb;
    std::vector<double> exps1 = bf1.GetExponents();
    std::vector<double> exps2 = bf2.GetExponents();
    std::vector<double> coefs1 = bf1.GetCoefficients();
    std::vector<double> coefs2 = bf2.GetCoefficients();
    for ( unsigned int m = 0; m < bf1.GetNumPrimitives(); m++ )
    {
      for ( unsigned int n = 0; n < bf2.GetNumPrimitives(); n++ )
      {
        mu = exps1[m] * exps2[n] / ( exps1[m] + exps2[n] );
        p = exps1[m] + exps2[n];
        xp = ( exps1[m] * atom1.GetX() + exps2[n] * atom2.GetX() ) / p;
        yp = ( exps1[m] * atom1.GetY() + exps2[n] * atom2.GetY() ) / p;
        zp = ( exps1[m] * atom1.GetZ() + exps2[n] * atom2.GetZ() ) / p;
        xpa = xp - atom1.GetX();
        ypa = yp - atom1.GetY();
        zpa = zp - atom1.GetZ();
        xpb = xp - atom2.GetX();
        ypb = yp - atom2.GetY();
        zpb = zp - atom2.GetZ();
        sum += coefs1[m]*coefs2[n]*bf1.GetNorm(m)*bf2.GetNorm(n)*
               Overlap_PRIM(mu,p,dx,dy,dz,i1,j1,k1,i2,j2,k2,xpa,ypa,zpa,xpb,ypb,zpb);
      }
    }
    return sum;
  }

  void Overlap ( const BUEHT::Atom & atom1, const BUEHT::Atom & atom2,
                    const BUEHT::BasisSet & basis1,
                    const BUEHT::BasisSet & basis2,
                    double* overlap_matrix )
  {
    int index = 0;
    double sum = 0.0;
    double dx(atom2.GetX()-atom1.GetX());
    double dy(atom2.GetY()-atom1.GetY());
    double dz(atom2.GetZ()-atom1.GetZ());
    for( unsigned int i = 0; i < basis1.GetNumShells(); i++ )
    {
      for ( int k = -1*basis1.GetL(i); k <= basis1.GetL(i); k++ ) // -m .. m
      {
        for ( unsigned int j = 0; j < basis2.GetNumShells(); j++ )
        {
          for ( int l = -1*basis2.GetL(j); l <= basis2.GetL(j); l++ ) // -m .. m
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
                                           basis1.GetBasisFunction(i),
                                           basis2.GetBasisFunction(j),
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
            overlap_matrix[index] = sum;
            index++;
          }
        }
      }
    }
  }
} // End of Namespace BUEHT

#endif
