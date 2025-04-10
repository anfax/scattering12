FC = ifort 
CMPLFLGS = -O3 -cm -w90 -w95 -ftz -static -ip 

LAPACK = -L/central/intel/mkl/lib/em64t -L/usr/lib64 -lguide -lmkl_solver -lmkl_lapack -lmkl_em64t -lpthread -Vaxlib
ROOT_DIR=/fourbody1/greene/tonzani/HHG_CPC/make/../src
LIB_DIR= /fourbody1/greene/tonzani/HHG_CPC/make/../src/../lib
ARPACK=/usr/local/lib/libarpack.a
SUPERLU_LIB = 
ARCR = ar -cr
EIGEN=eigen_pardiso.f90
COULOMB=call_coul.sub.f90
CERNLIB=-L/usr/local/cernlib/lib/ -lpacklib -lphtools -lmathlib -lkernlib
PARDISOSUB=pardiso.sub.f90
SPARSE_SUB=sparse_zerobased.f90
.SUFFIXES:
.SUFFIXES: .f .f90 .o
