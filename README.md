# pypercarlo
HyPerCarlo is going pythonic!

Idea here is to use the ALPS package http://alps.comp-phys.org/mediawiki/index.php/Main_Page and build application around that. No need to write everything from scratch as I tried in HyPerCarlo. HyPerCarlo will also be updated to use the ALPS C++ libraries.

## Structure of program
The program currently only has two classes, 

1. Spin - which allows the instantiation of vectors orientated in 3D space with an associated magnitude and typical vector operations (currently only the scalar product).

2. Simulation - which is the main controlling object and allows the user to implement a model of their choice and run Monte Carlo sims on it, this uses the ALPS library.

## TODO

1. Write methods to quickly visualise the spin lattice, need to check it is ordering properly in the more complicated models.

2. Test the different flip functions to make sure they are working correctly.

3. Write an Ewald summation module or method to call in dipole-dipole models.

## Planned simulations:

A variety ....
