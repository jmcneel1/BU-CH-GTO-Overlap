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

/*
  This is a struct to store angular momentum i,j,k alongside
  a normalization factor when defining the real spherical
  harmonics below
*/

struct RealSphericalHarmonic 
{
  double factor;
  int x, y, z;
}

/*
  Here we store the expansion of real spherical harmonics
  in terms of cartesians

  We order as m = -l .. 0 .. +l

  Also the order above for RealSphericalHarmonic holds ... x, y, z
*/

std::vector<RealSphericalHarmonic> RealSphericalHarmonics [] =
{
  {0.282094791773878,0,0,0}, // S
  {0.488602511902920,0,1,0}, // PY
  {0.488602511902920,0,0,1}, // PZ
  {0.488602511902920,1,0,0}, // PX
  {1.09254843059208,1,1,0},  // XY
  {1.09254843059208,0,1,1},  // YZ
  {0.630783130505040,0,0,2},  // Z2
  {-0.315391565252520,2,0,0},  // Z2
  {-0.315391565252520,0,2,0},  // Z2
  {1.09254843059208,1,0,1},  // XZ
  {0.546274215296040,2,0,0},  // X2-Y2
  {-0.546274215296040,0,2,0},  // X2-Y2
  {1.77013076977993,2,1,0},  // Y(3X2-Y2)
  {-0.590043589926644,0,3,0},  // Y(3X2-Y2)
  {2.89061144264055,1,1,1},  // xyz
  {1.82818319785786,0,1,2},  // YZ2
  {-0.457045799464466,2,1,0},  // YZ2
  {-0.457045799464466,0,3,0},  // YZ2
  {0.746352665180231,0,0,3},  // Z3
  {-1.11952899777035,2,0,1},  // Z3
  {-1.11952899777035,0,2,1},  // Z3
  {1.82818319785786,1,0,2},  // XZ2
  {-0.457045799464466,1,2,0},  // XZ2
  {-0.457045799464466,3,0,0},  // XZ2
  {1.44530572132028,2,0,1},  // Z(X2-Y2)
  {-1.44530572132028,0,2,1},  // Z(X2-Y2)
  {-1.77013076977993,1,2,0},  // X(X2-3Y2)
  {0.590043589926644,3,0,0},  // X(X2-3Y2)
  {2.50334294179670,3,1,0},  // XY(X2-Y2)
  {-2.50334294179670,1,3,0},  // XY(X2-Y2)
  {5.31039230933979,2,1,1},  // YZ(3*X2-Y2)
  {-1.77013076977993,0,3,1},  // YZ(3*X2-Y2)
  {5.67704817454536,1,1,2},  // XY(6*Z2-X2-Y2)
  {-0.946174695757560,3,1,0},  // XY(6*Z2-X2-Y2)
  {-0.946174695757560,1,3,0},  // XY(6*Z2-X2-Y2)
  {2.67618617422916,0,1,3},  // YZ(4Z2-3X2-3Y2)
  {-2.00713963067187,2,1,1},  // YZ(4Z2-3X2-3Y2)
  {-2.00713963067187,0,3,1},  // YZ(4Z2-3X2-3Y2)
  {0.84628437532163,0,0,4},  // 8Z^4-24X2Z2-24Y2Z2+3X4+3Y4+6X2Y2
  {-2.53885312596490,2,0,2},  // 8Z^4-24X2Z2-24Y2Z2+3X4+3Y4+6X2Y2
  {-2.53885312596490,0,2,2},  // 8Z^4-24X2Z2-24Y2Z2+3X4+3Y4+6X2Y2
  {0.317356640745613,4,0,0},  // 8Z^4-24X2Z2-24Y2Z2+3X4+3Y4+6X2Y2
  {0.317356640745613,0,4,0},  // 8Z^4-24X2Z2-24Y2Z2+3X4+3Y4+6X2Y2
  {0.634713281491226,2,2,0},  // 8Z^4-24X2Z2-24Y2Z2+3X4+3Y4+6X2Y2
}

} // End of namespace BUEHT

#endif
