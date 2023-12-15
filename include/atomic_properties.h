#ifndef _bueht_atomic_properties
#define _bueht_atomic_properties

/*

  James McNeely

  This file contains some useful arrays that
  store information about the atoms modeled 
  with BUEHT

  The individual arrays are described below
*/

namespace BUEHT
{

// Masses taken from IUPAC (2021)
// https://iupac.qmul.ac.uk/AtWt/

float bueht_atomic_masses [] = {
  1.008, // H
  4.002602, // He
  6.94, // Li
  9.0121831, // Be
  10.81, // B
  12.011, // C
  14.007, // N
  15.999, // O
  18.998403163, // F
  20.1797, // Ne
  22.98976928, // Na
  24.305, // Mg
  26.9815384, // Al
  28.085, // Si
  30.973762, // P
  32.06, // S
  35.45, // Cl
  39.95, // Ar
  39.0983, // K
  40.078, // Ca
  44.955907, // Sc
  47.867, // Ti
  50.9415, // V
  51.996, // Cr
  54.938043, //Mn
  55.845, // Fe
  58.933194, // Co
  58.6934, // Ni
  63.546, // Cu
  65.38, // Zn
  69.723, // Ga
  72.63, // Ge
  74.92160, // As
  78.97, // Se
  79.904, // Br
  83.798, // Kr
  85.4678, // Rb
  87.62, // Sr
  88.905838, // Y
  91.224, // Zr
  92.90637, // Mo
  97., // Tc
  101.07, // Ru
  102.90549, // Rh
  106.42, // Pd
  107.8682, // Ag
  112.414, // Cd
  114.818, // In
  118.71, // Sn
  121.760, // Sb
  127.60, // Te
  126.90447, // I
  131.29, // Xe
  132.905452, // Cs
  137.33, // Ba
  138.9055, // La
  140.116, // Ce
  140.90766, // Pr
  144.242, // Nd
  145., // Pm
  150.36, // Sm
  151.964, // Eu
  157.25, // Gd
  158.92535, // Tb
  162.500, // Dy
  164.93033, // Ho
  167.259, // Er
  168.93422, // Tm
  173.045, // Yb
  174.9668, // Lu
  178.49, // Hf
  180.94788, // Ta
  183.84, // W
  186.207, // Re
  190.23, // Os
  192.217, // Ir
  195.08, // Pt
  196.96657, // Au
  200.592, // Hg
  204.38, // Th
  207.2, // Pb
  208.98040, // Bi
  209., // Po
  210., // At
  222., // Rn
  223., // Fr
  226., // Ra
  227., // Ac
  232.0377, // Th
  231.03588, // Ta
  238.02891, // U
  237., // Np
  244., // Pu
  243., // Am
  247., // Cm
  247., // Bk
  251., // Cf
  252., // Es
  257., // Fm
  258., // Md
  259., // No
  262., // Lr
  267., // Rf
  270., // Db
  269., // Sg
  270., // Bh
  270., // Hs
  278., // Mt
  281., // Ds
  281., // Rg
  285., // Cn
  286., // Nh
  289., // Fl
  289., // Mc
  293., // Lv
  293., // Ts
  294. // Og
};

std::string bueht_atomic_symbols[] = {
      "H",
      "He",
      "Li",
      "Be",
      "B",
      "C",
      "N",
      "O",
      "F",
      "Ne",
      "Na",
      "Mg",
      "Al",
      "Si",
      "P",
      "S",
      "Cl",
      "Ar",
      "K",
      "Ca",
      "Sc",
      "Ti",
      "V",
      "Cr",
      "Mn",
      "Fe",
      "Co",
      "Ni",
      "Cu",
      "Zn",
      "Ga",
      "Ge",
      "As",
      "Se",
      "Br",
      "Kr",
      "Rb",
      "Sr",
      "Y",
      "Zr",
      "Nb",
      "Mo",
      "Tc",
      "Ru",
      "Rh",
      "Pd",
      "Ag",
      "Cd",
      "In",
      "Sn",
      "Sb",
      "Te",
      "I",
      "Xe",
      "Cs",
      "Ba",
      "La",
      "Ce",
      "Pr",
      "Nd",
      "Pm",
      "Sm",
      "Eu",
      "Gd",
      "Tb",
      "Dy",
      "Ho",
      "Er",
      "Tm",
      "Yb",
      "Lu",
      "Hf",
      "Ta",
      "W",
      "Re",
      "Os",
      "Ir",
      "Pt",
      "Au",
      "Hg",
      "Tl",
      "Pb",
      "Bi",
      "Po",
      "At",
      "Rn",
      "Fr",
      "Ra",
      "Ac",
      "Th",
      "Pa",
      "U",
      "Np",
      "Pu",
      "Am",
      "Cm",
      "Bk",
      "Cf",
      "Es"
      "Fm",
      "Md",
      "No",
      "Lr",
      "Rf",
      "Db",
      "Sg",
      "Bh",
      "Hs",
      "Mt",
      "Ds",
      "Rg",
      "Cn",
      "Nh",
      "Fl",
      "Mc",
      "Lv",
      "Ts",
      "Og"
    };

    std::string bueht_atomic_names[] = {
      "HYDROGEN", // 1
      "HELIUM", // 2
      "LITHIUM", // 3
      "BERYLLIUM", // 4
      "BORON", // 5
      "CARBON", // 6
      "NITROGEN", // 7
      "OXYGEN", // 8
      "FLUORINE", // 9 
      "NEON", // 10
      "SODIUM", // 11
      "MAGNESIUM", // 12
      "ALUMINUM", // 13
      "SILICON", // 14
      "PHOSPHORUS", // 15
      "SULFUR", // 16
      "CHLORINE", // 17 
      "ARGON", // 18
      "POTASSIUM", // 19
      "CALCIUM", // 20
      "SCANDIUM", // 21
      "TITANIUM", // 22
      "VANADIUM", // 23
      "CHROMIUM", // 24
      "MANGANESE", // 25
      "IRON", // 26
      "COBALT", // 27
      "NICKEL", // 28
      "COPPER", // 29
      "ZINC", // 30
      "GALLIUM", // 31
      "GERMANIUM", // 32
      "ARSENIC", // 33
      "SELENIUM", // 34
      "BROMINE", // 35
      "KRYPTON", // 36
      "RUBIDIUM", // 37
      "STRONTIUM", // 38
      "YTTRIUM", // 39
      "ZIRCONIUM", // 40
      "NIOBIUM", // 41
      "MOLYBDENUM", // 42
      "TECHNETIUM", // 43
      "RUTHENIUM", // 44
      "RHODIUM", // 45
      "PALLADIUM", // 46
      "SILVER", // 47
      "CADMIUM", // 48
      "INDIUM", // 49
      "TIN", // 50
      "ANTIMONY", // 51
      "TELLURIUM", // 52
      "IODINE", // 53
      "XENON", // 54
      "CAESIUM", // 55
      "BARIUM", // 56
      "LANTHANUM", // 57
      "CERIUM", // 58
      "PRASEODYMIUM", // 59
      "NEODYMIUM", // 60
      "PROMETHIUM", // 61
      "SAMARIUM", // 62
      "EUROPIUM", // 63
      "GADOLINIUM", // 64
      "TERBIUM", // 65
      "DYSPROSIUM", // 66
      "HOLMIUM", // 67
      "ERBIUM", // 68
      "THULIUM", // 69
      "YTTERBIUM", // 70
      "LUTETIUM", // 71
      "HAFNIUM", // 72
      "TANTALUM", // 73
      "TUNGSTEN", // 74
      "RHENIUM", // 75
      "OSMIUM", // 76
      "IRIDIUM", // 77
      "PLATINUM", // 78
      "GOLD", // 79
      "MERCURY", // 80
      "THALLIUM", // 81
      "LEAD", // 82
      "BISMUTH", // 83
      "POLONIUM", // 84
      "ASTATINE", // 85
      "RADON", // 86
      "FRANCIUM", // 87
      "RADIUM", // 88
      "ACTINIUM", // 89
      "THORIUM", // 90
      "PROTACTINIUM", // 91
      "URANIUM", // 92
      "NEPTUNIUM", // 93
      "PLUTONIUM", // 94
      "AMERICIUM", // 95
      "CURIUM", // 96
      "BERKELIUM", // 97
      "CALIFORNIUM", // 98
      "EINSTEINIUM", // 99
      "FERMIUM", // 100
      "MENDELEVIUM", // 101
      "NOBELIUM", // 102
      "LAWRENCIUM", // 103
      "RUTHERFORDIUM", // 104
      "DUBNIUM", // 105
      "SEABORGIUM", // 106
      "BOHRIUM", // 107
      "HASSIUM", // 108
      "MEITNERIUM", // 109
      "DARMSTADTIUM", // 110
      "ROENTGENIUM", // 111
      "COPERNICIUM", // 112
      "NIHONIUM", // 113
      "FLEROVIUM", // 114
      "MOSCOVIUM", // 115
      "LIVERMORIUM", // 116
      "TENNESSINE", // 117
      "OGANESSON" // 118
    };

} // End of namespace BUEHT

#endif
