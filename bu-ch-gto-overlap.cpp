#include <iostream>
#include <iomanip>
#include "atomic_properties.h"
#include "overlap.h"
#include "bueht_util.h"

/*

  James McNeely

  This is a simple tool that will calculate the overlap 
  integral between different atoms. It is intended
  to be used as a pedagogical tool perhaps by providing 
  complementary information to a Walsh diagram...
  
  All available shells for a given basis set will be printed.

  Thus will be prompted for the following:

  1: The atomic number of atom 1.
  2: The coordinates (X,Y,Z-Angstroem) of atom 1.
  4. The same information for atom 2
  5. The Basis Set (options in share folder)

  A small table will be printed that will show the resulting overlaps.

*/

int main ( int argc, char* argv[] )
{
  // l1 = L for atom 1, l2 = L for atom 2
  // n1 is the n quantum number for atom 1
  // n2 is the n quantum number for atom 2

  // First we assume that the share directory is in the same
  // folder as the executable...

  std::string basis_loc = BUEHT::GetBasesLocation(argv[0]);

  int atomicnum1, atomicnum2;
 
  // The coordinates (angstrom)

  double x1, x2, y1, y2, z1, z2;

  // The name of the basis set and a placeholder variable for streaming

  std::string basis_name, line;
  
  // Now we start interacting with the user

  std::cout << "Welcome!\nThis is a simple utility to calculate the overlap\n";
  std::cout << "integral between atoms.\n";
  std::cout << "Please enter the atomic number of center 1: ";
  std::cin >> atomicnum1;
  std::cout << "X-Coordinate atom 1 (Angstroem): ";
  std::cin >> x1;
  std::cout << "Y-Coordinate atom 1 (Angstroem): ";
  std::cin >> y1;
  std::cout << "Z-Coordinate atom 1 (Angstroem): ";
  std::cin >> z1;

  std::cout << "Please enter the atomic number of center 2: ";
  std::cin >> atomicnum2;
  std::cout << "X-Coordinate atom 2 (Angstroem): ";
  std::cin >> x2;
  std::cout << "Y-Coordinate atom 2 (Angstroem): ";
  std::cin >> y2;
  std::cout << "Z-Coordinate atom 2 (Angstroem): ";
  std::cin >> z2;

  std::cout << "Please enter the basis set name: ";
  std::cin >> basis_name;

  // Declare the atoms

  BUEHT::Atom atom1(atomicnum1,x1,y1,z1);
  BUEHT::Atom atom2(atomicnum2,x2,y2,z2);

  // Bring the name to uppercase

  BUEHT::StringToUpper(basis_name);

  // And now make sure it is supported

  int index = 0;
  bool basis_good = false;
  while ( index < *(&BUEHT::bueht_supported_basis_sets+1) - BUEHT::bueht_supported_basis_sets )
  {
    if ( basis_name == BUEHT::bueht_supported_basis_sets[index] )
    {
      basis_good = true;
      break;
    }
    index++;
  }

  if (!basis_good)
  {
    std::cout << basis_name << " not supported! Exiting..." << std::endl;
    std::exit(1);
  }
  
  // Read in the basis set.
  // The BasisSet class will perform all necessary checks.

  BUEHT::BasisSet basis1(basis_loc+basis_name+".bs",atomicnum1);
  BUEHT::BasisSet basis2(basis_loc+basis_name+".bs",atomicnum2);

  double overlap_matrix[basis1.GetDimensions() * basis2.GetDimensions()];

  std::cout << basis1.GetDimensions() << "\n";

  BUEHT::Overlap(atom1,atom2,basis1,basis2,overlap_matrix);

  /*

  // Check to make sure that this BF exists for atom 1.
  // If not, prompt the user again...

  while ( true )
  {
    if ( (*(BUEHT::bueht_param_ptr[atomicnum1-1]+count)).an == atomicnum1 )
    {
      if ((*(BUEHT::bueht_param_ptr[atomicnum1-1]+count)).l == l1 )
      {
        n1 = (*(BUEHT::bueht_param_ptr[atomicnum1-1]+count)).n;
        break;
      }
    }
    else
    {
      std::cout << "Invalid L for atom 1 with atomic number " << atomicnum1 << "...\n";
      std::cout << "Please enter the angular moments (L=0,1,2) for atom 1: ";
      std::cin >> l1;
      count = -1;
    }
    count = count + 1;
  }
  
  count = 0;
  while ( true )
  {
    if ( (*(BUEHT::bueht_param_ptr[atomicnum2-1]+count)).an == atomicnum2 )
    {
      if ((*(BUEHT::bueht_param_ptr[atomicnum2-1]+count)).l == l2 )
      {
        n2 = (*(BUEHT::bueht_param_ptr[atomicnum2-1]+count)).n;
        break;
      }
    }
    else
    {
      std::cout << "Invalid L for atom 2 with atomic number " << atomicnum2 << "...\n";
      std::cout << "Please enter the angular moments (L=0,1,2) for atom 2: ";
      std::cin >> l2;
      count = -1;
    }
    count = count + 1;
  }
  // Now we actually calculate the overlap
  std::vector<double> coord1 {x1,y1,z1};
  std::vector<double> coord2 {x2,y2,z2};
  BUEHT::Atom atom1(atomicnum1,coord1);
  BUEHT::Atom atom2(atomicnum2,coord2);

  // Note basis function parameters (zeta) were taken from Greg Landrum's 
  // YAeHMOP (https://github.com/greglandrum/yaehmop.git)

  BUEHT::BasisFunction bf1 (n1,l1,atomicnum1);
  BUEHT::BasisFunction bf2 (n2,l2,atomicnum2);
  if ( l1 == 0 && l2 == 0 ) // Two s-shells
  {
    double overlap = BUEHT::OverlapSS(bf1,bf2,atom1,atom2);
    std::cout << "\n";
    std::cout << std::setw(14) << BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"s" << std::endl;;
    std::cout << std::setw(7) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"s";
    std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap << "\n\n";
  }
  else if ( (l1 == 0 && l2 == 1) || (l1 == 1 && l2 == 0) ) // One S and one P shell
  {
    std::vector<double> overlap = BUEHT::OverlapSP(bf1,bf2,atom1,atom2);
    if ( l1 == 0 )
    {
      std::cout << "\n";
      std::cout << std::setw(14) << BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"px";
      std::cout << std::setw(7) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"py";
      std::cout << std::setw(7) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"pz" << std::endl;;
      std::cout << std::setw(7) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"s";
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[0];
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[1];
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[2] << "\n\n";
    }
    else
    {
      std::cout << "\n";
      std::cout << std::setw(14) << BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"s\n";
      std::cout << std::setw(7) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"px";
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[0] << std::endl;
      std::cout << std::setw(7) <<  BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"py";
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[1] << std::endl;
      std::cout << std::setw(7) <<  BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"pz";
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[2] << "\n\n";
    }
  }
  else if ( l1 == 1 && l2 == 1 )
  {
      std::vector<double> overlap = BUEHT::OverlapPP(bf1,bf2,atom1,atom2);
      std::cout << "\n";
      std::cout << std::setw(14) << BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"px";
      std::cout << std::setw(7) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"py";
      std::cout << std::setw(7) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"pz" << std::endl;;
      std::cout << std::setw(7) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"px";
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[0];
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[1];
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[2] << "\n";
      std::cout << std::setw(7) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"py";
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[3];
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[4];
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[5] << "\n";
      std::cout << std::setw(7) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"pz";
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[6];
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[7];
      std::cout << std::setw(7) << std::fixed << std::setprecision(3) << overlap[8] << "\n\n";
  }
  else if ( (l1 == 0 && l2 == 2) || (l1 == 2 && l2 == 0) ) // One S and one D shell
  {
    std::vector<double> overlap = BUEHT::OverlapSD(bf1,bf2,atom1,atom2);
    if ( l1 == 0 )
    {
      std::cout << "\n";
      std::cout << std::setw(20) << BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dxy";
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dyz";
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dz2";
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dxz";
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dx2-y2" << std::endl;;
      std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"s";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[0];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[1];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[2];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[3];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[4] << "\n\n";
    }
    else
    {
      std::cout << "\n";
      std::cout << std::setw(20) << BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"s"<<std::endl;
      std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dxy";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[0] << std::endl;
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dyz";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[1] << std::endl;
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dz2";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[2] << std::endl;
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dxz";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[3] << std::endl;
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dx2-y2";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[4] << "\n\n";
    }
  }
  else if ( (l1 == 1 && l2 == 2) || (l1 == 2 && l2 == 1) ) // One P and one D shell
  {
    std::vector<double> overlap = BUEHT::OverlapPD(bf1,bf2,atom1,atom2);
    if ( l1 == 1 )
    {
      std::cout << "\n";
      std::cout << std::setw(20) << BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dxy";
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dyz";
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dz2";
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dxz";
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dx2-y2" << std::endl;;
      std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"px";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[0];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[1];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[2];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[3];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[4] << std::endl;
      std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"py";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[5];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[6];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[7];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[8];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[9] << std::endl;
      std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"pz";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[10];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[11];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[12];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[13];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[14] << "\n\n";
    }
    else 
    {
      std::cout << "\n";
      std::cout << std::setw(20) << BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"px";
      std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"py";
      std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"pz" << std::endl;
      std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dxy";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[0];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[5];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[10] << std::endl;
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dyz";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[1];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[6];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[11] << std::endl;
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dz2";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[2];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[7];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[12] << std::endl;
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dxz";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[3];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[8];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[13] << std::endl;
      std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dx2-y2";
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[4];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[9];
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[14] << "\n\n";
    }
  }
  else if ( (l1 == 2 && l2 == 2) ) // Two D shells
  {
    std::vector<double> overlap = BUEHT::OverlapDD(bf1,bf2,atom1,atom2);
    std::cout << "\n";
    std::cout << std::setw(20) << BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dxy";
    std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dyz";
    std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dz2";
    std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dxz";
    std::cout << std::setw(10) <<  BUEHT::bueht_atomic_symbols[atomicnum2-1]+"-"+std::to_string(n2)+"dx2-y2" << std::endl;;
    std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dxy";
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[0];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[1];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[2];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[3];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[4] << std::endl;
    std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dyz";
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[5];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[6];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[7];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[8];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[9] << std::endl;
    std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dz2";
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[10];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[11];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[12];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[13];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[14] << std::endl;
    std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dxz";
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[15];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[16];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[17];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[18];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[19] << std::endl;
    std::cout << std::setw(10) << BUEHT::bueht_atomic_symbols[atomicnum1-1]+"-"+std::to_string(n1)+"dx2-y2";
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[20];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[21];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[22];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[23];
    std::cout << std::setw(10) << std::fixed << std::setprecision(3) << overlap[24] << std::endl << std::endl;
  }
  */
}
