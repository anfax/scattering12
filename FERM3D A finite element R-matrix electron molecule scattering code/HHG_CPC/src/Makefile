include ./make.inc

#****************************************************************
#****************************************************************

# Sources
SRCF90 = modules_3d.f90 gau_2D_tab.f90 over_calc.f90 gau_tab.f90 V_setup.f90 function_integr_3D.f90 solve.f90 kmatrix.f90 grid.f90 mat_el.f90 diago_sparse.f90 pleg.f90 struct_mat.f90 potential.f90 $(COULOMB) $(PARDISOSUB) $(SPARSE_SUB) ortog.f90 Integral_coulomb.f90 integr_1D.var_points.f90 integr_coul_fixedpts.f90 separate_integral.f90 $(EIGEN) photoionization.f90

SRCF77 = besselnew.f rgg.f Ay.gensub.f

NORMAL= sparse_eigen.f90 



# Get objects
OBJECT = $(SRCF77:.f=.o) $(SRCF90:.f90=.o)


OBJNORMAL = $(NORMAL:.f90=.o)

MYLIB = $(LIB_DIR)/liblocal.a

LIBNORMAL = $(LIB_DIR)/libnormal.a


#SPIENV =  $(ROOT_DIR)/SuperLUCustom/sp_ienv.o
#SLU = $(SUPERLU_LIB)
LIB1 = $(LIB_DIR)/libprop1d.a
LIBSPLINES = $(LIB_DIR)/libsplines_pot_interp.a


default: $(MYLIB) FERM3D.x 

# Make libraries
$(MYLIB): $(OBJECT)	
	$(ARCR) $(MYLIB) $(OBJECT)

$(LIBNORMAL): $(OBJNORMAL)
	$(ARCR) $(LIBNORMAL) $(OBJNORMAL)

# Pardiso Version
FERM3D.x: $(MYLIB) $(LIBNORMAL) $(SLU) $(LIB1) main.o
	export LD_LIBRARY_PATH
	#export OMP_NUM_THREADS=$(thr)
	$(FC) $(LOADFLGS) -o ../bin/FERM3D.x main.o $(MYLIB) $(LIBNORMAL) $(SPIENV)  $(LIB1) $(SLU) $(LIBSPLINES) $(CERNLIB) $(ARPACK) $(LAPACK) $(LINK_F91) 

# Make clean
clean: 
	rm *.o *.x *.a *.mod fort.*


# Compiling stuff
.f90.o:
	$(FC) $(CMPLFLGS) -c $<

.f.o:
	$(FC) $(CMPLFLGS) -c $<
