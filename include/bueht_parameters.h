#ifndef _bueht_parameters
#define _bueht_parameters

#include <bueht_constants.h>

namespace BUEHT
{

/*
  James McNeely

  The following definitions store parameters used by the program
  including basis function information.
  
*/

/*
  A paramater object is a "row" of the parameters
  listed below
*/

std::string bueht_supported_basis_sets[] = {
  "6-31G*",
  "DEF2-SV(P)",
  "DEF2-SVP",
  "DEF2-TZVP",
  "STO-3G"
};

} // End of namespace BUEHT

#endif
