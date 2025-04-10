!!!*************************************************************
! 文件/File: strtfunc.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: strtfunc.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      THIS FILE INCLUDES SUBROUTINES DEVOTED TO CALCULATION     C
C      OF THE MATRIX ELEMENTS ON THE STRETCHING  BASE FUNCTIONS  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : CALCULATES MATRIX ELEMENTS OF THE DIPOLE MOMENT OF THE 
C             : TRIATOMIC MOLECULE USING A NUMERICAL INTAGRATION 
C             : ON THE STRETCHING WAVE FUNCTIONS 
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################
      SUBROUTINE StretchMatrixElement 
C
c ---- Quantum number parameters for stretching wave functions ----- 
      include 'dipolsys.h'
      include 'molcul.h'
      include 'mu.h'

      REAL*8 PrimMorseMatrElements(0:4,2,2,0:MaxV1V3,0:MaxV1V3)
      REAL*8 StretchFunc(MaxNumberV1V3,MaxNvibMax,3,2)

      REAL*8  StretchMatrixElements(15,2,2,2,UpperNvib,UpperNvib)
      COMMON /DIPSTRETCH/ StretchMatrixElements


C--------------------------------------------------------------------C
C--------------------------------------------------------------------C
C --------------------- START OF THE SUBROUTINE ---------------------C

      WRITE(6,3001)
3001  FORMAT(1H0,' ABTOC.INT.INF   STRETCHING MATRIX ELEMENTS '/)


C----      READ EIGENVECTORS COEFITIENTS OF THE STRETCHING SYMMETR. WAVE FUNCTIONS
c      IF (SYMM) THEN
         CALL ReadStretchFunction(StretchFunc)
c      ENDIF


C-----     CALCULATION OF THE MATRIX ELEMENTS ON THE PRIMITIVE ------
C-----     UNVibETRIZED MORSE  FUNCTIONS                      ------
      CALL Calc_PrimMorseMatrElements
     1                        (PrimMorseMatrElements)


C----      CALCULATE MATRIX ELEMENTS OF THE Y_1^i*Y_3^j ON THE SYMMETR. 
C----      WAVE FUNCTIONS  <N,Gamma| Y_1^i*Y_3^j |N`,Gamma`>
C----      !!!! RIGHT HAND SIDE FUNCTION IS "C" ELECTR. STATE !!!!
      CALL Calc_SymmMorseMatrElements(StretchFunc,
     1                                PrimMorseMatrElements)

      RETURN 
      END 

C########################      StretchMatrixElement   #########################
C###########################################################################


C##########################################################################
C#             :
C# DATE        : 05.06.2000
C# AUTHOR      : Yurchenko
C# UPDATES     :
C# LANGUAGE    : FORTRAN
C# PURPOSE     : CALCULATION OF THE MATRIX ELEMENTS ON THE PRIMITIVE
C#             : UNVibETRIZED MORSE  FUNCTIONS 
C# SUBPROGRAMS : SCALIN, FMORSE, SimpsonIntegral
C#      CALLED :
C#             :
C###########################################################################
      SUBROUTINE  Calc_PrimMorseMatrElements
     1                        (PrimMorseMatrElements)  
c
      include 'dipolsys.h'
      include 'mu.h'

      REAL*8 PrimMorseMatrElements(0:4,2,2,0:MaxV1V3,0:MaxV1V3)

c ---- Quantum number parameters for stretching wave functions ----- 
      REAL*8 R12term,Rinterval,RStep,SCALIN,FMORSE,RE(2),AA(2)
      include 'molcul.h'
C
      INTEGER ISurf,NBond,V0,Npts,Vi,Vj,MaxV1V3Max
      REAL*8 ALPHA,Scalefactor,RminValue,SimpsonIntegral
      REAL*8 FnctsMorseL(0:NumPtsBond,0:MaxV1V3),
     1       FnctsMorseR(0:NumPtsBond,0:MaxV1V3),
     2       Phi_V_Phi(0:NumPtsBond),Y(0:NumPtsBond)

C--------------------------------------------------------------------C
C--------------------------------------------------------------------C
C --------------------- START OF THE SUBROUTINE ---------------------C

      WRITE(6,3001)
3001  FORMAT(1H0,' ABTOC.INT.INF   PRIMARY MATRIX ELEMENTS '/)


C ----- CALCULATE MATRIX ELEMENTS < Morse (A Nbond  V_i) | Y^k | Morse (C Nbond  V_j) > 
C ----- AND                       < Morse (B Nbond  V_i) | Y^k | Morse (C Nbond  V_j) > 


      RE(1) = RE12
      RE(2) = RE32

      AA(1) = AA1
      AA(2) = AA3


C
C CALCULATE K VALUES FOR THE TWO MORSE OSCILLATORS; STORE IN RMK
C

C---          MORBID C STATE     --- 
      DO NBond = 1,2       ! --- 1-2, 1-3 bond 

C ---- INTAGRATION INTERVAL RMax - Rmin
      Rinterval = Rmax(Nbond) - Rmin(Nbond)
      RminValue = Rmin(Nbond)
ccccccccccccccccccccccccc????????????????????????????????
      RStep = Rinterval/(NumPtsBond-1)
ccccccccccccccccccccccccc????????????????????????????????

C---- CALCULATE POWERS OF MORSE OSCILATOR FOR EVERY POINTS IN THE INTERVAL ----
      DO Npts=0,NumPtsBond,1
         Y(Npts)= 1 - DEXP(-AA(NBond)*(Npts*Rstep+RminValue-RE(NBond)))
      ENDDO 

      DO V0 = 0, NV1V3Max(3)
C ----- SCALE FACTOR -----
        ALPHA=2.0D+00*(Kappa(NBond,3)-V0)-1.0D+00
        Scalefactor = SCALIN(AA(NBond),ALPHA,V0)

C---- CALCULATE MORSE FUNCTIONS FOR EVERY POINTS IN THE INTERVAL ----
      DO Npts=0,NumPtsBond,1
         R12term=Npts*Rstep+RminValue
         FnctsMorseR(Npts,V0)=
     1        FMORSE(V0,
     2               R12term,
     3               Kappa(NBOND,3),
     4               RE0(NBond,3),
     5               AA(NBond),
     6               Scalefactor)
      ENDDO 

      ENDDO      !----- V   -----

C---          RENNER A-B STATES     --- 

      DO ISurf = 1,2       ! --- A,B

      DO V0 = 0, NV1V3Max(ISurf)
C ----- SCALE FACTOR -----
        ALPHA=2.0D+00*(Kappa(NBond,ISurf)-V0)-1.0D+00
        Scalefactor = SCALIN(AA(NBond),ALPHA,V0)

C---- CALCULATE MORSE FUNCTIONS FOR EVERY POINTS IN THE INTERVAL ----
      DO Npts=0,NumPtsBond,1
         R12term=Npts*Rstep+RminValue
         FnctsMorseL(Npts,V0)=
     1        FMORSE(V0,
     2               R12term,
     3               Kappa(NBOND,ISURF),
     4               RE0(NBond,ISurf),
     5               AA(NBond),
     6               Scalefactor)
      ENDDO 

      ENDDO      !----- V0   -----



      DO Vi = 0, NV1V3Max(ISurf)  ! ---  A,B  
      DO Vj = 0, NV1V3Max(3)  ! ---  C

C0     -------  <i|  1  |j> --------- 
C0
C0  --- Determine integration function Phi_V_Phi ---
       DO Npts = 0,NumPtsBond
          Phi_V_Phi(Npts) = FnctsMorseL(Npts,Vi)*
     1                      FnctsMorseR(Npts,Vj)
       ENDDO 

C0
C0      NUMERICAL INTAGRATION WITH SIMPSON'S RULE WITH RESPECT TO R BOND
C0  
        PrimMorseMatrElements(0,Nbond,ISurf,Vi,Vj) =
     1                   SimpsonIntegral(NumPtsBond,Rinterval,Phi_V_Phi)



C1    -------  <i|  Y^1  |j> -------------------------------------------
C1
C1  --- Determine integration function Phi_V_Phi ---
       DO Npts = 0,NumPtsBond

          Phi_V_Phi(Npts) = FnctsMorseL(Npts,Vi)*
     1                      Y(Npts)*
     2                      FnctsMorseR(Npts,Vj)

       ENDDO 

C1
C1      NUMERICAL INTAGRATION WITH SIMPSON'S RULE WITH RESPECT TO R BOND 
C1  
        PrimMorseMatrElements(1,Nbond,ISurf,Vi,Vj) =
     1                   SimpsonIntegral(NumPtsBond,Rinterval,Phi_V_Phi)

C2    -------  <i|  Y^2  |j> -------------------------------------------
C2
C2  --- Determine integration function Phi_V_Phi ---
       DO Npts = 0,NumPtsBond

          Phi_V_Phi(Npts) = FnctsMorseL(Npts,Vi)*
     1                        Y(Npts)**2*
     2                      FnctsMorseR(Npts,Vj)

       ENDDO 

C2
C2      NUMERICAL INTAGRATION WITH SIMPSON'S RULE WITH RESPECT TO R BOND 
C2  
        PrimMorseMatrElements(2,Nbond,ISurf,Vi,Vj) =
     1                   SimpsonIntegral(NumPtsBond,Rinterval,Phi_V_Phi)


C3    -------  <i|  Y^3  |j> -------------------------------------------
C3
C3  --- Determine integration function Phi_V_Phi ---
       DO Npts = 0,NumPtsBond

          Phi_V_Phi(Npts) = FnctsMorseL(Npts,Vi)*
     1                        Y(Npts)**3*
     2                      FnctsMorseR(Npts,Vj)

       ENDDO 

C3
C3      NUMERICAL INTAGRATION WITH SIMPSON'S RULE WITH RESPECT TO R BOND 
C3  
        PrimMorseMatrElements(3,Nbond,ISurf,Vi,Vj) =
     1                   SimpsonIntegral(NumPtsBond,Rinterval,Phi_V_Phi)

C4    -------  <i|  Y^4  |j> -------------------------------------------
C4
C4  --- Determine integration function Phi_V_Phi ---
       DO Npts = 0,NumPtsBond

          Phi_V_Phi(Npts) = FnctsMorseL(Npts,Vi)*
     1                        Y(Npts)**4*
     2                      FnctsMorseR(Npts,Vj)

       ENDDO 

C4 
C4      NUMERICAL INTAGRATION WITH SIMPSON'S RULE WITH RESPECT TO R BOND 
C4  
        PrimMorseMatrElements(4,Nbond,ISurf,Vi,Vj) =
     1                   SimpsonIntegral(NumPtsBond,Rinterval,Phi_V_Phi)


      ENDDO      !----- Vj   -----
      ENDDO      !----- Vi   -----


      ENDDO      !----- Nbond   -----
      ENDDO      !----- ISurf   -----

      RETURN
      END

C##################      Calc_PrimMorseMatrElements    #####################
C###########################################################################

C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : READ EIGENVECTORS COEFITIENTS OF THE STRETCHING SYMMETR.
C             : WAVE FUNCTIONS NEEDED TOCALCULATE THE MATRIX ELEMENTS OF THE DIPOLE MOMENT 
C             : OF THE TRIATOMIC MOLECULE USING THE MORBID SCHEME.
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################
      SUBROUTINE ReadStretchFunction(StretchFunc)
C     IMPLICIT REAL*8 (A-H,O-Z)
C     IMPLICIT LOGICAL (A-Z)

      include 'mu.h'
      include 'dipolsys.h'

      REAL*8 StretchFunc(MaxNumberV1V3,MaxNvibMax,3,2)
C
      CHARACTER*200 LINE
      CHARACTER*80 FILNAM
      INTEGER ISurf,GammaVib,NVib0,N,V1,V3,i0,MaxN,Vterm
      INTEGER IW1(MaxNumberV1V3),IW3(MaxNumberV1V3),MaxV1,MaxV3 
      REAL*8 VTEMP,CoefMax

C--------------------------------------------------------------------C
C--------------------------------------------------------------------C
C --------------------- START OF THE SUBROUTINE ---------------------C



C     OPEN FILE WITH LOWER STATE STRATCHING FUNCTION
C
      FILNAM=DirInput(1:NumDirLetters)//FileStretchAB
      OPEN(UNIT=93,FILE=FILNAM)



C
C     READ A1-B1-ELECTRONIC STATE STRATCHING BASIS FUNCTION - Stretchunc
      DO ISurf = 1,2       ! --- A,B - electronic states
      READ(93,1000) LINE
      READ(93,*) Kappa(1,ISurf),Kappa(2,ISurf)
C
C     SKIP 3 LINES
C
      READ(93,1000) LINE
      READ(93,1000) LINE
      READ(93,1000) LINE

C
C     READ NEXT LABELING LINES
C
      i0=0
      DO N = 0,NV1V3Max(ISurf) ! --- N = 0..(V1+V3)_max
      DO V1=0,N                       ! --- V1
         i0=i0+1
         READ(93,1011) IW1(i0),IW3(i0)
C     CHECKING INPUT DATA ???????????????
      ENDDO                ! --- V1
      ENDDO                ! --- N = 0..(V1+V3)_max
      MaxN = i0 

      DO GammaVib=1,2                          ! -- A1,B2 
      DO NVib0 = 1,Nvib(ISurf,GammaVib)  ! -- NVib0  

C     ORDER OF INDEXES SHOULD BE 01 - 10 - 02 - 11 - 20 - 03 - 12 - 21 - 30 ANS SO ON
C
C     -- CHECKING INPUT DATA 
            READ(93,*) VTEMP
            IF (ABS(INT(Vtemp)-Vtemp).LT.0.1E-8.AND.Vtemp.NE.0)  THEN  
             WRITE(*,*) 'ERROR: INPUT DATA DO NOT',
     1                  ' CORRESPOND TO MAX PARAM-Rs ',
     2                  'IN THE SUBROUTINE ReadStretchFunction' 
            STOP
            ENDIF

      StretchFunc(1,NVib0,ISurf,GammaVib)=VTEMP

      MaxV1 = 0
      MaxV3 = 0

      CoefMax = DABS(VTEMP)

      DO i0 = 2,MaxN ! --- i0 = 1..MaxN
          Vterm = 1
          DO N = 1,NV1V3Max(ISurf) ! --- N = 0..(V1+V3)_max
          DO V1=0,N                       ! --- V1
             V3=N-V1
             Vterm = Vterm +1 
             IF (IW1(i0).EQ.V1.AND.IW3(i0).EQ.V3) THEN
                READ(93,*) VTEMP
                StretchFunc(Vterm,NVib0,ISurf,GammaVib) = VTEMP
                IF ( CoefMax.LE.ABS(VTEMP) .AND. 
     2               MOD(V3+GammaVib,2).EQ.1 ) THEN 
                   CoefMax = ABS(VTEMP)      
                   MaxV1 = V1
                   MaxV3 = V3
                ENDIF
             ENDIF
          ENDDO                ! --- V1
          ENDDO                ! --- N = 0..(V1+V3)_max
      ENDDO                ! --- i0 = 1..MaxN


      NotateStretch(ISurf,GammaVib,NVib0) = MaxV1*100+MaxV3

      ENDDO                ! --- NVib0  
      ENDDO                ! --- A1,B2 
      ENDDO                ! --- A,B - electronic states

      CLOSE(93)



C
C     OPEN FILE WITH UPPER STATE STRATCHING FUNCTION
C

      FILNAM=DirInput(1:NumDirLetters)//FileStretchC
      OPEN(UNIT=93,FILE=FILNAM)
C
C     READ C-ELECTRONIC STATE STRATCHING BASIS FUNCTION - Stretchunc
      ISurf = 3 

C
C     SKIP FIRST LINE
C
      READ(93,1000) LINE

      READ(93,*) Kappa(1,ISurf),Kappa(2,ISurf)

      READ(93,1000) LINE
      READ(93,1000) LINE
      READ(93,1000) LINE


C
C     READ NEXT LABELING LINES
C
      i0=0
      DO N = 0,NV1V3Max(ISurf) ! --- N = 0..(V1+V3)_max
      DO V1=0,N                       ! --- V1
         i0=i0+1
         READ(93,*) IW1(i0),IW3(i0)
	   continue 
C     CHECKING INPUT DATA ???????????????
      ENDDO                ! --- V1
      ENDDO                ! --- N = 0..(V1+V3)_max


      DO GammaVib=1,2                          ! -- A1,B2 
      DO NVib0 = 1,Nvib(ISurf,GammaVib)  ! -- NVib0  

C     ORDER OF INDEXES SHOULD BE 01 - 10 - 02 - 11 - 20 - 03 - 12 - 21 - 30 ANS SO ON
C
C     -- CHECKING INPUT DATA 
         READ(93,*) VTEMP
c            IF (ABS(INT(Vtemp)-Vtemp).LT.0.1E-8.AND.Vtemp.NE.0)  THEN  
c             WRITE(*,*) 'ERROR: INPUT DATA DO NOT',
c     1                  ' CORRESPOND TO MAX PARAM-Rs',
c     2                  'IN THE SUBROUTINE ReadStretchFunction' 
c            STOP
c            ENDIF
         StretchFunc(1,NVib0,ISurf,GammaVib)=VTEMP

      MaxV1 = 0
      MaxV3 = 0

      CoefMax = ABS(VTEMP)


      DO i0 = 2,MaxN ! --- i0 = 1..MaxN
          Vterm = 1
          DO N = 1,NV1V3Max(ISurf) ! --- N = 0..(V1+V3)_max
          DO V1=0,N                       ! --- V1
             V3=N-V1
             Vterm = Vterm +1 
             IF (IW1(i0).EQ.V1.AND.IW3(i0).EQ.V3) THEN
                READ(93,*) VTEMP
                StretchFunc(Vterm,NVib0,ISurf,GammaVib) = VTEMP
                IF ( CoefMax.LE.ABS(VTEMP) .AND. 
     2               MOD(V3+GammaVib,2).EQ.1 ) THEN 
                   CoefMax = ABS(VTEMP)      
                   MaxV1 = V1
                   MaxV3 = V3
                ENDIF
             ENDIF
          ENDDO                ! --- V1
          ENDDO                ! --- N = 0..(V1+V3)_max
      ENDDO                ! --- i0 = 1..MaxN

      NotateStretch(ISurf,GammaVib,NVib0) = MaxV1*100+MaxV3

      ENDDO                ! --- NVib0  
      ENDDO                ! --- A1,B2 


      CLOSE(93)


1000  FORMAT(A100)
1011  FORMAT(2I12)


      RETURN

      END

C#########################      ReadStretchunction    ######################
C###########################################################################


C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : CALCULATE MATRIX ELEMENTS OF THE Y_1^i*Y_3^j ON THE SYMMETR. 
C             : WAVE FUNCTIONS  <N,Gamma| Y_1^i*Y_3^j |N`,Gamma`>
C             : !!!! RIGHT SIDE FUNCTION IS "C" ELECTR. STATE !!!!
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################
      SUBROUTINE Calc_SymmMorseMatrElements(StretchFunc,
     1                                      PrimMorseMatrElements)

      include 'dipolsys.h'
      include 'molcul.h'
      include 'mu.h'

      REAL*8 ElementTemp

      REAL*8  StretchMatrixElements(15,2,2,2,UpperNvib,UpperNvib)
      COMMON /DIPSTRETCH/ StretchMatrixElements


      REAL*8 StretchFunc(MaxNumberV1V3,MaxNvibMax,3,2)
      REAL*8 PrimMorseMatrElements(0:4,2,2,0:MaxV1V3,0:MaxV1V3)

      INTEGER ISumPower,i0,ISurfL,GammaVibL,GammaVibR,
     1        NVibL,NVibR,NL,NR,VL1,VL3,VR1,VR3,
     2        IndexY1(15),IndexY3(15),ISurfR,iL,iR

C--------------------------------------------------------------------C
C--------------------------------------------------------------------C
C --------------------- START OF THE SUBROUTINE ---------------------C

      WRITE(6,3001)
3001  FORMAT(1H0,' ABTOC.INT.INF   FINAL MATRIX ELEMENTS '/)


C    THE ARRAYS IndexY CONTAIN THE COEFFICIENTS OF THE FOLLOWING
C    POWERS OF THE Y1 AND Y3 QUANTITIES
C
C    I  Used Symm    MUABY(I)
C  --------------------------------
C    1     1           1
C    2     2           Y1
C    3                 Y3
C    4     3           Y1**2
C    5                 Y3**2
C    6     4           Y1*Y3
C    7     5           Y1**3
C    8                 Y3**3
C    9     6           Y1**2*Y3
C   10                 Y1*Y3**2
C   11     7           Y1**4
C   12                 Y3**4
C   13     8           Y1**3*Y3
C   14                 Y1*Y3**3
C   15     9           Y1**2*Y3**2

      IndexY1( 1) = 0 
      IndexY1( 2) = 1
      IndexY1( 3) = 0
      IndexY1( 4) = 2
      IndexY1( 5) = 0
      IndexY1( 6) = 1
      IndexY1( 7) = 3
      IndexY1( 8) = 0
      IndexY1( 9) = 2
      IndexY1(10) = 1
      IndexY1(11) = 4
      IndexY1(12) = 0
      IndexY1(13) = 3
      IndexY1(14) = 1
      IndexY1(15) = 2

      IndexY3( 1) = 0 
      IndexY3( 2) = 0
      IndexY3( 3) = 1
      IndexY3( 4) = 0
      IndexY3( 5) = 2
      IndexY3( 6) = 1
      IndexY3( 7) = 0
      IndexY3( 8) = 3
      IndexY3( 9) = 1
      IndexY3(10) = 2
      IndexY3(11) = 0
      IndexY3(12) = 4
      IndexY3(13) = 1
      IndexY3(14) = 3
      IndexY3(15) = 2



C    
C     LOOP OVER DIFFERENT POWERS OF THE Y_1^i*Y_3^j
C
      DO i0 = 1,15             ! ---- i0 

C
C     LOOP OVER DIFFERENT QUANTUM NUMBERS OF THE SYMMTR. STRETCH. EIGEN. FUNCT. OF THE H_str 
C
      ISurfR = 3
      DO ISurfL = 1,2       ! --- A,B - electronic states  Left  < |
      DO GammaVibL=1,2                          ! -- A1,B2 Left  < |
      DO GammaVibR=1,2                          ! -- A1,B2 Right | >
      DO NVibL = 1,Nvib(ISurfL,GammaVibL)  ! -- NVibLeft < | 
      DO NVibR = 1,Nvib(3,GammaVibR)  ! -- NVibRight| >
C
C     LOOP OVER DIFFERENT QUANTUM NUMBERS OF THE PRIMM. STRETCH. EIGEN. FUNCT. OF THE H_morse
C
C     ORDER OF INDEXES SHOULD BE 01 - 10 - 02 - 11 - 20 - 03 - 12 - 21 - 30 ANS SO ON

      iL = 0
      ElementTemp = 0 
      DO NL = 0,NV1V3Max(ISurfL) ! --- N = 0..(V1+V3)_max Left 
      DO VL1= 0,NL               ! --- V1Left
         VL3 = NL-VL1            ! --- V3Left
         iL=iL+1
      iR = 0
      DO NR = 0,NV1V3Max(ISurfR) ! --- N = 0..(V1+V3)_max Right
      DO VR1= 0,NR               ! --- V1Right
         VR3 = NR-VR1            ! --- V3Right
         iR=iR+1
         ElementTemp = ElementTemp + 
     1   StretchFunc(iL,NVibL,ISurfL,GammaVibL)*
     2   PrimMorseMatrElements(IndexY1(i0),1,ISurfL,VL1,VR1)*  
     3   PrimMorseMatrElements(IndexY3(i0),2,ISurfL,VL3,VR3)*
     4   StretchFunc(iR,NVibR,ISurfR,GammaVibR)

	   continue 


      ENDDO                ! --- V1Right
      ENDDO                ! --- N = 0..(V1+V3)_max Right
      ENDDO                ! --- V1Left
      ENDDO                ! --- N = 0..(V1+V3)_max Left 
      StretchMatrixElements
     1            (i0,ISurfL,GammaVibL,GammaVibR,NVibL,NVibR) = 
     1            ElementTemp
         continue  


      ENDDO                ! --- ! -- NVibRight| >
      ENDDO                ! --- ! -- NVibLeft < | 
      ENDDO                ! --- ! -- A1,B2 Right | >
      ENDDO                ! --- ! -- A1,B2 Left  < |

      ENDDO                ! --- A,B - electronic states  Left  < |
      ENDDO                ! ---- i0 - powers Y_1^i*Y_3^j


      RETURN



      END
C##################      Calc_SymmMorseMatrElements   ######################
C###########################################################################



