!!!*************************************************************
! 文件/File: pqx2xyz.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: pqx2xyz.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE PQX2XYZ
      IMPLICIT REAL*8 (A-H,O-Z)

C             :
C DATE        : 05.06.2000
C AUTHOR      : 
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : SETS UP THE RHO DEPENDENT FUNCTIONS NEEDED TO
C             : CALCULATE THE ROTATION-VIBRATION ENERGIES OF A
C             : TRIATOMIC MOLECULE USING THE MORBID SCHEME.
C SUBPROGRAMS :
C      CALLED :
C             :
C

C
C
C
      REAL*8 EPSC1,EPSC2,CORHOE,HSTEP

      INTEGER I,L,LMax

      REAL*8 DMMACY(15),
     1       DMMACZ(15),
     2       DMMBCX(15),SQ2

      include 'dipolsys.h'
      include 'value.h'
      include 'molcul.h'
      include 'bcoeff.h'
      include 'crcoef.h'
c      include 'dimen.h'
      include 'mu.h'


      REAL*8 CoefMuRhoY(15,MaxNumRhoPoints),
     1       CoefMuRhoZ(15,MaxNumRhoPoints),
     1       CoefMuRhoX(15,MaxNumRhoPoints)

      CHARACTER*80 FILNAM

C
C---------------------------------------------------C
C---------------------------------------------------C
C ------------ START OF THE SUBROUTINE -------------C

      M = M1 + M2 + M3

      U1=M1*(M3+M2)*RE12*RE12
      U3=M3*(M1+M2)*RE32*RE32
      U13=M1*M3*RE12*RE32
      SQ2 = DSQRT(2.0D+00)
      V=M2*(M1+M2+M3)*RE12*RE32
      EPSC1=(U1-U3)/DSQRT((U1+U3)**2-4.0E0*U13**2)
      EPSC2=DSQRT((U1+U3-2.0E0*U13)/(U1+U3+2.0E0*U13))

      HSTEP =  RHOMAX/NumPtsRho
C
C     START LOOP OVER RHO
C
      DO 1000 I=1,NumPtsRho
      RHO=DFLOAT(I)*HSTEP
C
C
      EPS=0.5D0*RHO+EPSC1*ATAN(EPSC2*TAN(0.5D0*RHO))
      CR=COS(RHO)
      SR=SIN(RHO)
      CSE=COS(EPS)
      SNE=SIN(EPS)
      CRE=CR*CSE+SR*SNE
      SRE=SR*CSE-CR*SNE
      EPSP=(U1+U13*CR)/(U1+U3+2.0E0*U13*CR)
      EPSPP=U13*SR*(U1-U3)/(U1+U3+2.0E0*U13*CR)**2
      EPSPPP=U13*(U1-U3)*(CR/(U1+U3+2.0E0*U13*CR)**2
     1     +4.0E0*U13*SR*SR/(U1+U3+2.0E0*U13*CR)**3)

C
C     GENERATE THE COEFFICIENTS IN THE EXPANSIONS OF THE S I'S
C     IN TERMS OF THE DELTA R I'S.
C
      CALL BCOEF
C
C     GENERATE THE COEFFICIENTS IN THE EXPANSIONS OF COS RHO BAR
C     IN TERMS OF THE DELTA R I'S.
C
      CALL COSRB
C
C     GENERATE THE COEFFICIENTS IN THE EXPANSIONS OF THE MU
C     TENSOR ELEMENTS AND THE VIBRATIONAL ANGULAR MOMENTA IN
C     TERMS OF THE DELTA R I'S.
C
      CALL MUTENS ( DMMACY , DMMACZ , DMMBCX )
C


      DO L = 1,15
      CoefMuRhoX(L,I)=DMMBCX(L)
      CoefMuRhoY(L,I)=DMMACY(L)
      CoefMuRhoZ(L,I)=DMMACZ(L)
      ENDDO 

C
1000  CONTINUE
C

      OPEN(UNIT=NFileMu,STATUS='SCRATCH',FORM='UNFORMATTED')

      DO L = 1,15
        WRITE(NFileMu) ( CoefMuRhoX(L,I), I = 1,NumPtsRho )
        WRITE(NFileMu) ( CoefMuRhoY(L,I), I = 1,NumPtsRho )
        WRITE(NFileMu) ( CoefMuRhoZ(L,I), I = 1,NumPtsRho )
      ENDDO 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC??????????????????????????
c      FILNAM=DirInput(1:NumDirLetters)//'mu_x.dat'
c      OPEN(UNIT=23,FILE=FILNAM)
c
c      DO I = 1,NumPtsRho
c        WRITE(23,3111) I,( CoefMuRhoX(L,I), L = 1,15 )
c      ENDDO 
c
c3111  FORMAT(I5,15(D16.8))
c      CLOSE(23)
c
c      FILNAM=DirInput(1:NumDirLetters)//'mu_y.dat'
c      OPEN(UNIT=23,FILE=FILNAM)
c
c      DO I = 1,NumPtsRho
c        WRITE(23,3111) I,( CoefMuRhoy(L,I), L = 1,15 )
c      ENDDO 
c
c      CLOSE(23)
c
c      FILNAM=DirInput(1:NumDirLetters)//'mu_z.dat'
c      OPEN(UNIT=23,FILE=FILNAM)
c
c      DO I = 1,NumPtsRho
c        WRITE(23,3111) I,( CoefMuRhoz(L,I), L = 1,15 )
c      ENDDO 
c
c      CLOSE(23)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC??????????????????????????

      REWIND ( NFileMu )

      WRITE(6,3001)
3001  FORMAT(1H0,' ABTOC.INT.INF   SET UP RHO-DEPENDED DIPOLE FUNCT.'/)

      RETURN

3101  FORMAT('1',5X,' Sym.   parameters ')
3102  FORMAT('1',5X,' UnSym. parameters ')
3103  FORMAT('1',5X,F12.8,27(D16.8))
3104  FORMAT('1',5X,F12.8,18(D16.8))




      END

