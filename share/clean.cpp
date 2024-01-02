#include <sstream>
#include <iostream>
#include <fstream>
#include "atomic_properties.h"

int main ()
{
  std::ifstream inFile("SARC2-DKH-QZVP.bs");
  std::ofstream outFile("SARC2-DKH-QZVP.tmp.bs");
  std::string line, tmp, element;
  std::stringstream ss;
  int index;
  while ( getline(inFile,line) )
  {
    if ( line[0] == 'N' )
    {
      ss.clear(); ss.str("");
      ss << line;
      ss >> tmp >> element;
      for ( unsigned int i = 0; i < 118; i++ )
      {
        if ( BUEHT::bueht_atomic_symbols[i] == element ) 
        {
          index = i;
          break;
        }
      }
      outFile << BUEHT::bueht_atomic_names[index] << "\n";
      std::cout << element << "\n";
    }
    else
    {
      if ( line.find("END\\") != std::string::npos )
      {
        outFile << "END\n";
      }
      else
      {
        outFile << line << "\n";
      }
    }
  }
  inFile.close();
  outFile.close();
}
