#ifndef _bueht_basis_function
#define _bueht_basis_function

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <bueht_parameters.h>
#include <math.h>
#include <bueht_util.h>
#include <vector>

/*
  James McNeely

  This BasisFunction class stores relevant information for a shell.

  Private Variables include:

  2. The angular momentum qN (myL)
  3. The contraction coefficients. These are for a normalized spherical BFs. 
     ie... c1^2 <x1|x1> + c2^2 <x2|x2> + c3^2 <x3|x3> + 2c1c3<x1|x3> +
           2c1c2 <x2|x1> + 2c2c3 <x3|x2> = 1 for three primitives...
  4. The Zetas for each primitive
  5. The number of primitives.

*/

namespace BUEHT
{

class BasisFunction
{
  public:

    BasisFunction () = default;

    BasisFunction ( int l, int num_primitive,
                   const std::vector<double> & coeffs, 
                   const std::vector<double> & zeta )
    {
      myL = l;
      myNorms.resize(coeffs.size());
      // First check to make sure that the arguments have the correct structure
      if ( (num_primitive == coeffs.size()) && (num_primitive == zeta.size()) )
      {
        this->myNumPrimitives = num_primitive;
        this->myZetas = zeta;
        this->myCoefficients = coeffs;
        for ( unsigned int i = 0; i < coeffs.size(); i++ )
        {
          // N = (2*(2*zeta)^(3/4)/pi^(1/4))*sqrt(2^l/(2l+1)!!)*sqrt(2*zeta)^l
          myNorms[i] = std::pow(bueht_pi,-0.25)*
                       std::pow(2.e0,1.75+(double)l)*
                       std::pow(BUEHT::DoubleFactorial(2*l+1),-0.5)*
                       std::pow(myZetas[i],0.75+(double)l/2.e0);
        }
      }
      else
      {
        std::cout << "Incompatible primitize count and vector size! Exiting...\n";
        std::exit(1);
      }
    }

    double GetExponent ( const int & index ) const
    {
      if  ( (index >= 0) && (index < myNumPrimitives) )
      {
        return myZetas[index];
      }
      else
      {
        std::cout << "Invalid index for AO coefficient vector! Exiting...\n";
        std::exit(1);
      }
    }

    double GetCoefficient ( const int & index ) const
    {
      if  ( (index >= 0) && (index < myNumPrimitives) )
      {
        return myCoefficients[index];
      }
      else
      {
        std::cout << "Invalid index for AO coefficient vector! Exiting...\n";
        std::exit(1);
      }
    }

    std::vector<double> GetCoefficients () const
    {
      return myCoefficients;
    }

    std::vector<double> GetExponents () const
    {
      return myZetas;
    }

    int GetL () const { return myL; }

    unsigned int GetNumPrimitives () const { return myNumPrimitives; }

    double GetNorm ( const int & index ) const 
    { 
      if  ( (index >= 0) && (index < myNumPrimitives) )
      {
        return myNorms[index];
      }
      else
      {
        std::cout << "Invalid index for AO Norm vector! Exiting...\n";
        std::exit(1);
      }
    }

    std::vector<double> GetNorms () const { return myNorms; }


    /*
      This overload is currently not used, but will be worked on in the 
      future, so it is kept for the time being.

      James McNeely - 12/15/23
    */

    friend std::ostream & operator << ( std::ostream & buffer, 
                                   const BUEHT::BasisFunction & basis_function )
    {
      buffer << "L: " << basis_function.myL << "\n";
      buffer << "Number of Primitives: " << basis_function.myNumPrimitives << "\n";
      for ( unsigned int i = 0; i < basis_function.myNumPrimitives; i++ )
      {
        buffer << "{" << basis_function.myCoefficients[i] << ", ";
        buffer << basis_function.myZetas[i] << "}";
        if ( i == ( basis_function.myNumPrimitives - 1 ) )
        {
          buffer << "\n";
        }
        else
        {
          buffer << ", ";
        }
      }
      return buffer;
    }

  private:
    int myL;
    unsigned int myNumPrimitives;
    std::vector<double> myCoefficients;
    std::vector<double> myZetas;
    std::vector<double> myNorms;

};

}

#endif
