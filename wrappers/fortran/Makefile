all:
	g++ -c ccode.cpp -Irichdem/include
	gfortran -c fortran.f90
	gfortran fortran.o ccode.o -lstdc++
