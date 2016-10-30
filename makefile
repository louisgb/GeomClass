#-- compiler and math libraries
f90c = ifort
libflags = -mkl
# f90c = gfortran
# libflags = -llapack

#-- specify path to the math libraries if needed
# libpath = -L/usr/lib/lapack

abmat.exe : abmat.o frompolyrate.o
	$(f90c) $(libpath) $(libflags) -o abmat.exe abmat.o frompolyrate.o
abmat.o : abmat.f90
	$(f90c) -c abmat.f90
frompolyrate.o : frommstor/frompolyrate.f90
	$(f90c) -c frommstor/frompolyrate.f90
clean:
	rm abmat.o frompolyrate.o
