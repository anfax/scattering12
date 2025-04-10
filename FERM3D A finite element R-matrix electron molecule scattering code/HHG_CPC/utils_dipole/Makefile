include ./make.inc

#****************************************************************
#dewfault :: normal version
#drago_iter :: version that uses iterative solver (beta version!)
#****************************************************************

SRCF90 = modules_3d.f90 gau3d.f90 lin_3drmatsparse.f90 pleg.f90 adda.f90 call_coul.sub.f90 gau_2D_tab.f90 Integral_coulomb.f90 integr_1D.var_points.f90 integr_coul_fixedpts.f90  R_matrix_post_production.f90 multipl_norm.f90 kmatrix.perturbed.f90 Hydrogen_perturbed.f90 integr_1D.f90 call_bes.sub.f90

SRCF77 = besselnew.f Ay.gensub.f


OBJECT = $(SRCF77:.f=.o) $(SRCF90:.f90=.o)



MYLIB = $(LIB_DIR)/liblocal_pp.a



#LIBSPLINES = $(ROOT_DIR)/splines/libsplines_pot_interp.a
LIBSPLINES = $(LIB_DIR)/libsplines_pot_interp.a

#DBS = $(ROOT_DIR)/splines/dbs3ine90.o

default: $(MYLIB) R_mat_post_proc.dipole.x 

$(MYLIB): $(OBJECT)	
	ar cr $(MYLIB) $(OBJECT)

R_mat_post_proc.dipole.x: $(MYLIB)  main.o
	$(FC) $(LOADFLGS) -o ../bin/R_mat_post_proc.dipole.x main.o $(MYLIB) $(LIBSPLINES) $(CERNLIB) $(LAPACK)  

CM.x: prog_CM.f90
	 $(FC) $(LOADFLGS) -o ../bin/CM.x prog_CM.f90

clean: 
	rm *.o *.x *.a *.mod fort.*

.f90.o:
	$(FC) $(CMPLFLGS) -c $<

.f.o:
	$(FC) $(CMPLFLGS) -c $<
