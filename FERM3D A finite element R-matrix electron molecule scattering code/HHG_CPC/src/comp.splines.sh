#/usr/bin/csh
ifort -O3 -cm -w90 -w95 -ftz -static -ip  -c bspline90_22.f90
ifort -O3 -cm -w90 -w95 -ftz -static -ip  -c dbs3in_pot_interp.f90
ar -cr /fourbody1/greene/tonzani/HHG_CPC/lib/libsplines_pot_interp.a dbs3in_pot_interp.o bspline90_22.o
rm *.o
cp *.mod ../
