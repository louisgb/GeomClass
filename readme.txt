This package contains tools to displace Cartesian coordinates of a molecule according to changes in internal coordinates. The iterative algorithm is described in Pulay et al., J. Am. Chem. Soc. 1979, 101, 2550 (DOI: 10.1021/ja00504a009).

--- Usage ---

* Compile the Fortran code by executing
  $ make
  or 
  $ ifort -mkl -o abmat.exe abmat.f90 frommstor/frompolyrate.f90
* Prepare the following input files and put them in a directory:
  cart0.txt: starting Cartesian coordinates to be displaced
  rdef0.txt: definition of redundant internal coordinates
  nrdef0.txt: definition of nonredundant internal coordinates, as linear combinations of redundant ones (or identical to redundant ones)
  nrdispl0.txt: displacement of nonredundant internal coordinates
* Execute the code by running
  $ python [directory_containing_code/]displace_int2cart.py [directory_containing_input_files] >out.txt
  where the parts in [] are optional. If directory_containing_input_files is not given, the code will search the current directory. 
* Check the output file (out.txt in the last example)

--- System requirements and dependencies ---

Linux (tested on CentOS release 6.8)
Intel Fortran Compiler with MKL libraries (tested with ifort version 13.1.3)
Python 2.7 (tested with version 2.7.8)
NumPy (tested with version 1.11.1)

--- Components ---

[displace_int2cart.py](python)
Displace the Cartesian coordinates according to changes in internal coordinates.
Need to call abmat.exe (see below)
- INPUT -
(In all the input files, empty lines or lines starting with a hash (#) are ignored. See the example inputs for details of format.)
cart0.txt : cartesian coordinates of molecule with mass. The first line is "[number of atoms] [irun]" where if irun=0 then abmat.exe will be called and the A and B matrices will be generated but the iterations to determine Cartesian displacements will not run. Each of the following lines is "[atom symbol] [mass] [x] [y] [z]". 
rdef0.txt : Definition of redundant internals.
nrdef0.txt : Definition of nonredundant internals in terms of redundant ones.
nrdispl0.txt : Displacement of internal coordinates. Each line is "[NR internal index] [displacement]"
- INTERMEDIATE OUTPUT -
See INPUT of [abmat.exe] below. 
- OUTPUT -
(standard output) : updated Cartesian coordinates in the format of "[atom symbol] [mass] [x] [y] [z]". 


[abmat.exe](f90; compile by "make" or "ifort -mkl -o abmat.exe abmat.f90 frommstor/frompolyrate.f90")
Calculate Wilson B and A matrices. (Better be called by displace_int2cart.py with irun=0 (see above).)
- INPUT -
cart.txt : cartesian coordinates of molecule with mass. The first line is the number of atoms. Each of the following lines is "[atom symbol] [mass] [x] [y] [z]". 
ibl.txt : Definition of bond lengths as outputted by displace_int2cart.py
iba.txt : Definition of bond angles. Format similar to above
ito.txt : Definition of torsions. Format similar to above
iob.txt : Definition of oop bends. Format similar to above
iod.txt : Definition of oop distances. Format similar to above
intnrdef.txt : Definition of nonredundant internals in terms of redundant internals defined by the above i*.txt
- OUTPUT -
bmatr.txt : B matrix of redundant internals. The first line is the two dimensions of the matrix (NIntR * 3N). Following is the matrix. 
bmat.txt : B matrix of nonredundant internals. The first line is the two dimensions of the matrix (3N-6 * 3N). Following is the matrix. 
amat.txt : A matrix of nonredundant internals. The first line is the two dimensions of the matrix (3N * 3N-6). Following is the matrix. 

--- Notes ---
---- Updated Oct. 2015:
----- Added support for new internal coordinate type: out-of-plane distance.
----- Currently support 5 types of coordinates: bond length, bend, torsion, oop bend, oop distance. 
----------------------------
