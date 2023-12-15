#ifndef _bueht_basis_set
#define _bueht_basis_set

#include "basis_function.h"
#include "atomic_properties.h"
#include "bueht_parameters.h"
#include <sstream>
#include <map>

/*

  James McNeely

  This is the BasisSet class.

  This is simply a vector of BasisFunctions, and contains
  some routines to get dimensions of the basis set.

  In addition, functions are provided to get information
  about the individual basis functions themselves.

*/

namespace BUEHT
{

class BasisSet
{

  public:

    /*
      This constructor reads the basis set for an atom from a file in the share folder.
    */

    BasisSet( const std::string & fname, int atomicnum )
    {
      // Variable declarations
      char l; 
      int l_int, tmp;
      int nprim;
      double exp, coeff;
      std::string line;
      std::string element = BUEHT::bueht_atomic_names[atomicnum-1];

      // Open the buffer and declare a stringstream buffer for later
      std::ifstream inFile(fname);
      std::stringstream sstream;
      if ( !inFile.good() )
      {
        std::cout << "Couldn't open basis set file! Exiting...\n";
        std::exit(1);
      }

      // Go to the element passed in...
      while ( !inFile.eof() )
      {
        getline(inFile,line);
        if ( line.find(element) != std::string::npos )
        {
          break;
        }
      }

      if ( !inFile.eof() )
      {
        getline(inFile,line);
        // In the basis set files in share, all atoms have "END"
        // after the BS info is complete...
        while ( line.find("END") == std::string::npos )
        {
          sstream.clear(); sstream.str(std::string());
          sstream << line;
          sstream >> l >> nprim;

          l_int = BUEHT::bueht_angular_momentum[l];

          std::vector<double> exps(nprim);
          std::vector<double> coefs(nprim);

          for ( unsigned int i = 0; i < nprim; i++ )
          {
            getline(inFile,line);
            sstream.clear(); sstream.str(std::string());
            sstream << line;
            sstream >> tmp >> exps[i] >> coefs[i];
          }

          BUEHT::BasisFunction bf(l_int,nprim,coefs,exps);
          myBasisSet.push_back(bf);
          getline(inFile,line);
        }
      }
      else
      {
        std::cout << "Invalid Basis for " << element << "! Exiting...\n";
        std::exit(1);
      }
      inFile.close();
    }

    int GetNumShells () const
    {
      return myBasisSet.size();
    }

    int GetDimensions () const
    {
      int dimensions(0);
      for ( unsigned int i = 0; i < myBasisSet.size(); i++ )
      {
        dimensions += ( 2 * myBasisSet[i].GetL() + 1 );
      }
      return dimensions;
    }

    int GetL ( const int & index ) const
    {
      if ( index < myBasisSet.size() )
      {
        return myBasisSet[index].GetL();
      }
      else
      {
        std::cout << "Invalid Basis Set Index in GetL! Exiting...\n";
        std::exit(1);
      }
    }

    // Return the canonical angular momentum char ...
    // Useful for printing.

    char GetL_Char ( const int & index ) const
    {
      int l = myBasisSet[index].GetL();
      if ( l == 0 ) return 'S';
      else if ( l == 1 ) return 'P';
      else if ( l == 2 ) return 'D';
      else if ( l == 3 ) return 'F';
      else if ( l == 4 ) return 'G';
      else if ( l == 5 ) return 'H';
      else return 'I';
    }

    std::vector<double> GetCoefficients ( int index ) const
    {
      if ( index < myBasisSet.size() )
      {
        return myBasisSet[index].GetCoefficients();
      }
      else
      {
        std::cout << "Invalid Basis Set Index in GetCoefficients! Exiting...\n";
        std::exit(1);
      }
    }

    std::vector<double> GetExponents ( int index ) const
    {
      if ( index < myBasisSet.size() )
      {
        return myBasisSet[index].GetExponents();
      }
      else
      {
        std::cout << "Invalid Basis Set Index in GetExponents! Exiting...\n";
        std::exit(1);
      }
    } 

    std::vector<double> GetNorms ( int index ) const
    {
      if ( index < myBasisSet.size() )
      {
        return myBasisSet[index].GetNorms();
      }
      else
      {
        std::cout << "Invalid Basis Set Index in GetNorms! Exiting...\n";
        std::exit(1);
      }
    }

    BUEHT::BasisFunction GetBasisFunction ( int index ) const
    {
      if ( index < myBasisSet.size() )
      {
        return myBasisSet[index];
      }
      else
      {
        std::cout << "Invalid Basis Set Index in GetBasisFunction! Exiting...\n";
        std::exit(1);
      }
    }

  private:
    std::vector<BUEHT::BasisFunction> myBasisSet;
};

}

#endif
