#/usr/bin/csh
make clean
ifort -O2 -c modules_3d.f90
ifort -O2 -c gau3d.f90
ifort -O2 -c Hydrogen_perturbed.f90

make 
make CM.x
make clean
