#ifndef _bueht_util
#define _bueht_util

#include "bueht_constants.h"
#include "bueht_parameters.h"

/*
  James McNeely

  Some small utility functions used throughout the program.

  So far, this included a function to calculate the factorial
  and a binomial coefficient.
  
*/

namespace BUEHT
{

  int CartesianExpansionLength (int l, int m)
  {
    int index = l*l + l - m;
    return (BUEHT::RealSphericalHarmonicsPtr[index+1] - BUEHT::RealSphericalHarmonicsPtr[index]);
  }

  int factorial ( int n )
  {
    int result = 1;
    for ( int j = 2; j <= n; j++ )
    {
      result = result * j;
    }
    return result;
  }

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
