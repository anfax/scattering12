subroutines = subroutines/bkmp2.f subroutines/bcpsi.f90 subroutines/angle.f90 subroutines/legendre.f90 \
              subroutines/smallF.f90 subroutines/propagation.f90 subroutines/damping.f90 subroutines/momentGrid.f90 \
              subroutines/kinetic.f90
modules = modules/constants.f90
main: $(subroutines) main.f95
	gfortran -m64 main.f95 $(subroutines) $(modules) -o main.exe -L/opt/acml5.3.1/gfortran64/lib -lacml

clean:
	rm *.out
