abmat.exe : abmat.o frompolyrate.o
	ifort -mkl -o abmat.exe abmat.o frompolyrate.o
abmat.o : abmat.f90
	ifort -mkl -c abmat.f90
frompolyrate.o : frommstor/frompolyrate.f90
	ifort -mkl -c frommstor/frompolyrate.f90
clean:
	rm abmat.o frompolyrate.o
