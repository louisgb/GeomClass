This package contains tools to displace Cartesian coordinates of a molecule according to changes in internal coordinates. The iterative algorithm is described in Pulay et al., J. Am. Chem. Soc. 1979, 101, 2550 (DOI: 10.1021/ja00504a009).

Currently support 5 types of coordinates: 
* bond length, specified by 2 atom indexes i-j
* bend, specified by 3 atom indexes i-j-k
* torsion, specified by 4 atom indexes, i-j-k-l. 
* out-of-plane bend, specified by 4 atom indexes, i-j-k-l. Defined as angle between bond i-l and plane j-l-k.
* out-of-plane distance, specified by 4 atom indexes, i-j-k-l. Defined as distance of central atom l from plane i-j-k.

--- System requirements and dependencies ---

* Linux (tested on CentOS release 6.8)
* Fortran compiler with linear algebra libraries (tested with ifort 13.1.3 with MKL and gfortran 4.8.0 with LAPACK)
* Python 2.7 (tested with version 2.7.8)
* NumPy (tested with version 1.11.1)

--- Usage ---

* Modify makefile to specify a Fortran compiler, if neccesary
* Compile the Fortran code by executing
  $ make
  (Dollar sign represents the terminal prompt.)
* Prepare the following input files and put them in a directory:
  cart0.txt: starting Cartesian coordinates to be displaced
  rdef0.txt: definition of redundant internal coordinates
  nrdef0.txt: definition of nonredundant internal coordinates, as linear combinations of redundant ones (or identical to redundant ones)
  nrdispl0.txt: displacement of nonredundant internal coordinates
* Execute the code by running
  $ python [directory_containing_code/]main.py [directory_containing_input_files] >out.txt
  where the parts in [] are optional. If directory_containing_input_files is not given, the code will search the current directory for input files. 
* Check the output file (out.txt in the last example)

For advanced users:
Write your own wrapper to use the GeomClass class. (Read Documentation.txt and use main.py as an example.)

--- Components ---

[main.py](Python)
Wrapper to displace Cartesian coordinates according to changes in internal coordinates.
- INPUT -
(In all the input files, empty lines and lines starting with a hash (#) are ignored. See the example inputs for details of format.)
  cart0.txt : Cartesian coordinates of molecule with atomic masses. 
  rdef0.txt : Definition of redundant internals.
  nrdef0.txt : Definition of nonredundant internals in terms of redundant ones.
  nrdispl0.txt : Displacement of internal coordinates. 
- INTERMEDIATE OUTPUT -
  See INPUT of [abmat.exe] below. 
- OUTPUT -
  (standard output) : updated Cartesian coordinates and convergence info. 


[geomclass.py](Python)
A class for molecular geometry. See Documentation.txt or run "$ pydoc geomclass" for more details. 
(Need to call abmat.exe [see below])


[geomclass_aux.py](Python)
Some auxillary functions.


[abmat.exe](Fortran 90)
Calculate Wilson B and A matrices. 
- INPUT -
  cart.txt : cartesian coordinates of molecule with mass. The first line is the number of atoms. Each of the following lines is "[atom symbol] [mass] [x] [y] [z]". 
  ibl.txt : Definition of bond lengths. The first line is the number of bond lengths. The following lines are pairs of atom indices of each bond. 
  iba.txt : Definition of bond angles. Format similar to above.
  ito.txt : Definition of torsions. Format similar to above.
  iob.txt : Definition of oop bends. Format similar to above.
  iod.txt : Definition of oop distances. Format similar to above.
  intnrdef.txt : Definition of nonredundant internals in terms of redundant internals defined by the above i*.txt
- OUTPUT -
  bmatr.txt : B matrix of redundant internals. The first line is the two dimensions of the matrix (NIntR * 3N). Following is the matrix. 
  bmat.txt : B matrix of nonredundant internals. The first line is the two dimensions of the matrix (3N-6 * 3N). Following is the matrix. 
  amat.txt : A matrix of nonredundant internals. The first line is the two dimensions of the matrix (3N * 3N-6). Following is the matrix. 


