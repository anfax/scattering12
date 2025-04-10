#
# NOTE: Only on some RedHat Linux systems the Absoft f77 compiler needs
#       the special option '-s' to compile the f77 version of the  demo
#       codes!! All others work fine.
# 
#------------------------------------------------------------------------------
#
# flags for FORTRAN 77 compiler
#

f77comp   = g77
f77flags  = 

#
#------------------------------------------------------------------------------
#
# flags for FORTRAN90 compiler
#

f90comp   = ifc
f90flags  = 



# *****************************************************************************

allf90targets =	dbsnake90 dbsinte90 dbsdere90 dbs1gde90 dbs2ine90 dbs2dre90 \
                dbs2gde90 dbs3ine90 dbs3dre90 dbs3gde90

allf90: 	$(allf90targets)

# *****************************************************************************

allf77targets = dbsnake dbsinte dbsdere dbs1gde dbs2ine dbs2dre dbs2gde \
                dbs3ine dbs3dre dbs3gde

allf77: 	$(allf77targets)

# *****************************************************************************

all: 	allf77 allf90

# *****************************************************************************

clean: 	
	rm -f *.o *.mod *~ *.a

# *****************************************************************************

distclean:	clean 	
	rm -f $(allf77targets) $(allf90targets)

# *****************************************************************************



# *****************************************************************************

dbsnake = dbsnake.o dlamch.o bspline_22.o

dbsnake: $(dbsnake)
	$(f77comp) $(f77flags) $(dbsnake) -o dbsnake

# *****************************************************************************

dbsnake90 = bspline90_22.o dbsnake90.o 

dbsnake90: $(dbsnake90)
	$(f90comp) $(f90flags) $(dbsnake90) -o dbsnake90

# *****************************************************************************

dbsinte = dbsinte.o dlamch.o bspline_22.o

dbsinte: $(dbsinte)
	$(f77comp) $(f77flags) $(dbsinte) -o dbsinte

# *****************************************************************************

dbsinte90 = bspline90_22.o dbsinte90.o 

dbsinte90: $(dbsinte90)
	$(f90comp) $(f90flags) $(dbsinte90) -o dbsinte90

# *****************************************************************************

dbsdere = dbsdere.o dlamch.o bspline_22.o

dbsdere: $(dbsdere)
	$(f77comp) $(f77flags) $(dbsdere) -o dbsdere

# *****************************************************************************

dbsdere90 = bspline90_22.o dbsdere90.o 

dbsdere90: $(dbsdere90)
	$(f90comp) $(f90flags) $(dbsdere90) -o dbsdere90

# *****************************************************************************

dbs1gde = dbs1gde.o dlamch.o bspline_22.o

dbs1gde: $(dbs1gde)
	$(f77comp) $(f77flags) $(dbs1gde) -o dbs1gde

# *****************************************************************************

dbs1gde90 = bspline90_22.o dbs1gde90.o 

dbs1gde90: $(dbs1gde90)
	$(f90comp) $(f90flags) $(dbs1gde90) -o dbs1gde90

# *****************************************************************************

dbs2ine = dbs2ine.o dlamch.o bspline_22.o

dbs2ine: $(dbs2ine)
	$(f77comp) $(f77flags) $(dbs2ine) -o dbs2ine

# *****************************************************************************

dbs2ine90 = bspline90_22.o dbs2ine90.o 

dbs2ine90: $(dbs2ine90)
	$(f90comp) $(f90flags) $(dbs2ine90) -o dbs2ine90

# *****************************************************************************

dbs2dre = dbs2dre.o dlamch.o bspline_22.o

dbs2dre: $(dbs2dre)
	$(f77comp) $(f77flags) $(dbs2dre) -o dbs2dre

# *****************************************************************************

dbs2dre90 = bspline90_22.o dbs2dre90.o 

dbs2dre90: $(dbs2dre90)
	$(f90comp) $(f90flags) $(dbs2dre90) -o dbs2dre90

# *****************************************************************************

dbs2gde = dbs2gde.o dlamch.o bspline_22.o

dbs2gde: $(dbs2gde)
	$(f77comp) $(f77flags) $(dbs2gde) -o dbs2gde

# *****************************************************************************

dbs2gde90 = bspline90_22.o dbs2gde90.o 

dbs2gde90: $(dbs2gde90)
	$(f90comp) $(f90flags) $(dbs2gde90) -o dbs2gde90

# *****************************************************************************

dbs3ine = dbs3ine.o dlamch.o bspline_22.o

dbs3ine: $(dbs3ine)
	$(f77comp) $(f77flags) $(dbs3ine) -o dbs3ine

# *****************************************************************************

dbs3ine90 = bspline90_22.o dbs3ine90.o 

dbs3ine90: $(dbs3ine90)
	$(f90comp) $(f90flags) $(dbs3ine90) -o dbs3ine90

# *****************************************************************************

dbs3dre = dbs3dre.o dlamch.o bspline_22.o

dbs3dre: $(dbs3dre)
	$(f77comp) $(f77flags) $(dbs3dre) -o dbs3dre

# *****************************************************************************

dbs3dre90 = bspline90_22.o dbs3dre90.o 

dbs3dre90: $(dbs3dre90)
	$(f90comp) $(f90flags) $(dbs3dre90) -o dbs3dre90

# *****************************************************************************

dbs3gde = dbs3gde.o dlamch.o bspline_22.o

dbs3gde: $(dbs3gde)
	$(f77comp) $(f77flags) $(dbs3gde) -o dbs3gde

# *****************************************************************************

dbs3gde90 = bspline90_22.o dbs3gde90.o 

dbs3gde90: $(dbs3gde90)
	$(f90comp) $(f90flags) $(dbs3gde90) -o dbs3gde90

# *****************************************************************************


.SUFFIXES: .o .f .f90

.f.o : ; $(f77comp) $(f77flags) -c $*.f
.f90.o : ; $(f90comp) $(f90flags) -c $*.f90

