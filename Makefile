HOMEPATH=$(PWD)
MODLOC=/usr2/postdoc/piskulic/privatemodules

check_dist:
	mkdir -p bin
	gfortran -fopenmp -O3 fortran/check_dist.f90 -o bin/check_dist.exe
