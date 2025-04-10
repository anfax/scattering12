#/usr/bin/csh
ifc -c eigensparse.f90
ifc -c sparsemultiply.f90
ifc -c eigenstuff.f90
ar cr libprop1d.a eigensparse.o sparsemultiply.o eigenstuff.o
