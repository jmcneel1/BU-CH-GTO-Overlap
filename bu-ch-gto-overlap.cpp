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

  // First we assume that the share directory is in the same
  // folder as the executable...

  std::string basis_loc = BUEHT::GetBasesLocation(argv[0]);

  // The atomic numbers of the two atomic centers

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

  BUEHT::Atom atom1(atomicnum1,x1,y1,z1,'A');
  BUEHT::Atom atom2(atomicnum2,x2,y2,z2,'A');

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

  BUEHT::Overlap(atom1,atom2,basis1,basis2,overlap_matrix);

  /*
    Now we print out the results
  */

   int shell_num = 0, m_index = 0, m_index2 = 0;

  // First print out the first column

  std::cout << std::setw(5) << ' ';

  for ( unsigned int i = 0; i < basis1.GetDimensions(); i++ )
  {
    std::cout << std::setw(7) << basis1.GetL_Char(shell_num);
    if ( ( m_index + 1 ) == ( 2 * basis1.GetL(shell_num) + 1 ) )
    {
      shell_num++;
      m_index = 0;
    }
    else
    {
      m_index++;
    }
  }

  std::cout << "\n";
  shell_num = 0;
  m_index = 0;

  for ( unsigned int i = 0; i < basis2.GetDimensions(); i++ )
  {
    std::cout << std::setw(5) << basis2.GetL_Char(shell_num);
    if ( ( m_index + 1 ) == ( 2 * basis1.GetL(shell_num) + 1 ) )
    {
      shell_num++;
      m_index = 0;
    }
    else
    {
      m_index++;
    }
    for ( unsigned int j = 0; j < basis1.GetDimensions(); j++ )
    {
      std::cout << std::setw(7) << std::fixed << std::setprecision(3)
                << overlap_matrix[i*basis1.GetDimensions()+j];
    }
    std::cout << "\n";
  }

  return 0; // HOOZAH!!

}
