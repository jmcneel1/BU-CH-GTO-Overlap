# BU-CH-GTO-Overlap
Simple Program to display overlap integrals between GTO shells on atom pairs.

Handles S-I shells, although benchmarks were only performed to H-Shells...

All modern C++ compilers should successfully compile the program. Because the program uses header files, no Makefile is provided here.

Simply run (in the case of GNU g++)

g++ -c -Iinclude/ bu-ch-gto-overlap.cpp

g++ -o bu-ch-gto-overlap bu-ch-gto-overlap.o

This makes an executable called bu-ch-gto-overlap

To run, simply use:

./bu-ch-gto-overlap
