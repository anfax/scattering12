#/usr/bin/csh
cd /fourbody1/greene/tonzani/HHG_CPC/make/../src
make clean
cd splines/
make clean
chmod +x comp.splines.sh
./comp.splines.sh
cp *.mod ../
cd ../
#cd SuperLUCustom
#chmod +x comp.SLU.sh
#./comp.SLU.sh
#cd ../
cd eigen_sparse
make clean
chmod +x comp.eigen.sh
./comp.eigen.sh
cd ../
ifort -O3 -cm -w90 -w95 -ftz -static -ip  -c modules_3d.f90
gmake
make clean
