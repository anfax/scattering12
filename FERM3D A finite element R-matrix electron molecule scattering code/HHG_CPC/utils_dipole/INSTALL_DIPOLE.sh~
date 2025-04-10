#/usr/bin/csh
cd /fourbody1/greene/tonzani/HHG_CPC/make/../src/../utils_dipole
make clean
ifort -O3 -cm -w90 -w95 -ftz -static -ip  -c modules_3d.f90
ifort -O3 -cm -w90 -w95 -ftz -static -ip  -c gau3d.f90
ifort -O3 -cm -w90 -w95 -ftz -static -ip  -c gau_2d_tab.f90
ifort -O3 -cm -w90 -w95 -ftz -static -ip  -c Hydrogen_perturbed.f90
make 
make CM.x
make clean
