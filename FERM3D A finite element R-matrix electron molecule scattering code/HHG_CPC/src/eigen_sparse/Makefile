include ../make.inc

SUPERLU = $(HOME)/SuperLU/superlu_alpha.a

SOURCE = eigenstuff.f90 eigenfew.f90 efield.f90 coefmask.f90 ham2d.f90 \
eigensparse.f90 sparsemultiply.f90 eigen2d.f90 chebysparse.f90 \
eigensparse1.f90 chebysparse1.f90 makebessel.f90
OBJECT = $(SOURCE:.f90=.o)
PROPLIB1 = libprop1d.a

default: $(PROPLIB1) testsparse prop1d prop2d test2d cheby2 cheby1 \
colin fd2 postfd2

$(PROPLIB1): $(OBJECT)
	ar cr $(PROPLIB1) $(OBJECT)
	ranlib $(PROPLIB1)

colin: $(PROPLIB1) colin.o
	$(FC) $(LOADFLGS) -o xco colin.o $(FEB_LIB) $(PROPLIB1) $(ARPACK) $(LAPACK) $(FEB_LIB) $(SUPERLU)

prop1d: $(PROPLIB1) prop1d.o
	$(FC) $(LOADFLGS) -o prop1d prop1d.o $(FEB_LIB) $(PROPLIB1) $(ARPACK) $(LAPACK)

prop2d: $(PROPLIB1) prop2d.o
	$(FC) $(LOADFLGS) -o prop2d prop2d.o $(FEB_LIB) $(PROPLIB1) $(ARPACK) $(LAPACK) $(FEB_LIB) $(SUPERLU)

cheby2: $(PROPLIB1) cheby2.o
	$(FC) $(LOADFLGS) -o xcheb2 cheby2.o $(FEB_LIB) $(PROPLIB1) $(ARPACK) $(LAPACK) $(FEB_LIB) $(SUPERLU)

cheby1: $(PROPLIB1) cheby1.o
	$(FC) $(LOADFLGS) -o xcheb1 cheby1.o $(FEB_LIB) $(PROPLIB1) $(ARPACK) $(LAPACK) $(FEB_LIB) $(SUPERLU)

test2d: $(PROPLIB1) test2d.o
	$(FC) $(LOADFLGS) -o test2d test2d.o $(FEB_LIB) $(PROPLIB1) $(ARPACK) $(LAPACK)

testsparse: $(PROPLIB1) testsparse.o
	$(FC) $(LOADFLGS) -o testsparse testsparse.o $(FEB_LIB) $(PROPLIB1) $(ARPACK) $(SUPERLU) $(LAPACK)

dndrv4: dndrv4.o 
	$(FC) $(LOADFLGS) -o dndrv4 dndrv4.o $(ARPACK) $(LAPACK)

fd2: $(PROPLIB1) fd2.o
	$(FC) $(LOADFLGS) -o fd2 fd2.o $(PROPLIB1) $(ARPACK)

postfd2: postfd2.o
	$(FC) $(LOADFLGS) -o postfd2 postfd2.o


.f90.o: $(FEB_MOD)
	$(FC) $(CMPLFLGS) -c $< $(VERBOSE)

clean:
	rm *.o
