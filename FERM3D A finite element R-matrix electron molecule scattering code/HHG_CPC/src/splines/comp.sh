#/usr/bin/csh
ifc -c bspline90_22.f90
ifc -c dbs3in_pot_interp.f90
ar cr libsplines_pot_interp.a dbs3in_pot_interp.o bspline90_22.o
cp *.mod ../
