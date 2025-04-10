# Makefile used for AMD OPTERON COMPUTING SYSTEMS, 
# Model TYAN GX28 (B2882).
# Operating system: SuSE LINUX PROFESSIONAL 9. 
# The compiler is pgi: the Portland Group 6.2 Compilers suite.

# MAKEFILE:
    FC = pgf90
    FP =  pgf77
    MPIF90 = mpif90
    CPP=cpp

    LIB    = -L/opt/acml/pgi64/lib -lacml -lacml_mv

    DEFINES= -DPARALLEL

# All ARCHs for debugging/profiling thus overriding defaults:
#
#   MYFLGS = -g
#   MYFLGS = -O2 -pg
#
    FFLAGS = $(MYFLGS) 

.SUFFIXES: .o .f .f90 

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90
.f.o:
	$(FP) $(FFLAGS) -c $*.f
clean:
	/bin/rm -f *.o *.bak core *.mod *~

realclean:
	/bin/rm -f *.o *.bak core *.exe *.mod *~

OBJ =   common_blocks.o close_coupling.o logd.o \
	vpot.o vfitt.o lib.o logdvp.o \
	var_ph.o derivs.o derivs1.o daprop_lr.o 

SRC =   common_blocks.f90 close_coupling.f90 logd.f90 \
	vpot.f90 vfitt.f90 lib.f logdvp.f90 \
	var_ph.f90 derivs.f90 derivs1.f90 daprop_lr.f90 
# 
serial: serial.f90 $(OBJ) $(SRC)
	$(FC) $(FFLAGS) $(LIB) -o cc.exe $(OBJ) serial.f90

parallel: parallel.f90 $(OBJ) $(SRC)
	$(MPIF90) $(FFLAGS) $(LIB) -o cc_mpi.exe parallel.f90 $(OBJ) 
#
serial.f90: main.cpp.f90
	$(CPP) -P -C main.cpp.f90  > serial.f90

parallel.f90: main.cpp.f90
	$(CPP) $(DEFINES) -P -C main.cpp.f90 > parallel.f90

