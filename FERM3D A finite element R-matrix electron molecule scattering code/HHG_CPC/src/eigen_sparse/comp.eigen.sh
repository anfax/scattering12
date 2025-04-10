#/usr/bin/csh
ifort -O3 -cm -w90 -w95 -ftz -static -ip  -c eigensparse.f90
ifort -O3 -cm -w90 -w95 -ftz -static -ip  -c sparsemultiply.f90
ifort -O3 -cm -w90 -w95 -ftz -static -ip  -c eigenstuff.f90
ar -cr /fourbody1/greene/tonzani/HHG_CPC/make/../src/../lib/libprop1d.a eigensparse.o sparsemultiply.o eigenstuff.o
rm *.o
