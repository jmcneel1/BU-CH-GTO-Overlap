#ifndef _bueht_basis_set
#define _bueht_basis_set

#include "basis_function.h"
#include "atomic_properties.h"
#include "bueht_parameters.h"
#include <sstream>
#include <map>

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
      char l;
      int l_int, tmp;
      int nprim;
      double exp, coeff;
      std::string line;
      std::string element = BUEHT::bueht_atomic_names[atomicnum-1];
      std::ifstream inFile(fname);
      std::stringstream sstream;
      if ( !inFile.good() )
      {
        std::cout << "Couldn't open basis set file! Exiting...\n";
        std::exit(1);
      }
      std::cout << element << "\n";
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
      return myBasisSet[index].GetL();
    }

    std::vector<double> GetCoefficients ( int index ) const
    {
      return myBasisSet[index].GetCoefficients();
    }

    std::vector<double> GetExponents ( int index ) const
    {
      return myBasisSet[index].GetExponents();
    } 

    std::vector<double> GetNorms ( int index ) const
    {
      return myBasisSet[index].GetNorms();
    }

  private:
    std::vector<BUEHT::BasisFunction> myBasisSet;
};

}

#endif
