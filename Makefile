HOMEPATH=$(PWD)

check_dist:
	mkdir -p bin
	gfortran -fopenmp -O3 fortran/check_distances.f90 -o bin/check_dist.exe
