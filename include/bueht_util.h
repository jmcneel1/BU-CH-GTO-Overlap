#ifndef _bueht_util
#define _bueht_util

#include "bueht_constants.h"
#include "bueht_parameters.h"

/*
  James McNeely

  Some small utility functions used throughout the program.
  
*/

namespace BUEHT
{

  /*

    This function returns the number of terms in the Cartesion expansion of a 
    given Real Spherical Harmonic.

    The angular momentum (L) and azimuthal number (M) are required.

    The function makes use of the "size" member of the RealSphericalHarmonic
    struct ( in bueht_parameters.h )

  */

  int CartesianExpansionLength (int l, int m)
  {
    int index = l*l + l - m;
    return BUEHT::RealSphericalHarmonics[BUEHT::RealSphericalHarmonicsPtr[index]].length;
  }

  /*
    n! = n*(n-1)*(n-2)...2*1
    Obviously, in this case 1 is returned if n < 2
  */

  int factorial ( int n )
  {
    int result = 1;
    for ( int j = 2; j <= n; j++ )
    {
      result = result * j;
    }
    return result;
  }

  /*
    ( n )       n!
    |   | =  --------
    ( k )    k!(n-k)!
  */

  double binomial_coeff ( int n, int k )
  {
    if ( ( 0 <= k ) && ( k < n ) )
    {
      return factorial(n)/(factorial(k)*factorial(n-k));
    }
    else
    {
      return 0;
    }
  }

  /*
    Brings a given string to capital letters.

    We node that the function toupper(char c) leaves 
    character unchanced unless c is a lowercase letter.
  */

  void StringToUpper ( std::string & caps )
  {
    int i = 0;
    char c;
    while (caps[i])
    {
      c = caps[i];
      caps[i] = toupper(c);
      i++;
    }
  }


  /*

    This is a little function to get the path of the 
    executable. The executable here will normally
    be argv[0].

    We use the compiler to determine if a Windows
    machine is being used, in which case the path separator
    much be altered.

  */

  std::string GetBasesLocation ( char* executable )
  {
    std::string str_executable(executable);
    #if defined(_WIN32) || defined(WIN32)
      int pos = str_executable.rfind('\\');
    #else
      int pos = str_executable.rfind('/');
    #endif

    std::string path = str_executable.substr(0,pos+1);

    return path+"share/";
  }

  /*

    (n)!! = 
      ODD N : n(n-2)(n-4)...(3)(1)
      EVEN N: n(n-2)(n-4)...(4)(2)
 
  */

  int DoubleFactorial ( const int & arg )
  {
    unsigned int total = 1;
    if ( (arg%2) == 0 )
    {
      for ( unsigned int i = 1; i <= arg/2; i++ ) total*=(2*i);
    }
    else
    {
      for ( unsigned int i = 1; i <= (arg+1)/2; i++ ) total*=(2*i-1);
    }
    return total;
  }
  
}

#endif
