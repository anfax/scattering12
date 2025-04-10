!!!*************************************************************
! 文件/File: abcrea.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: abcrea.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : 
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : READ DIPOLE MOMENT AB INICIO DATA 
C             : CALCULATED IN THE (P-Q-X) COORDINATE SISTEM 
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################

      SUBROUTINE ReadInputData
      LOGICAL POTSYM
      CHARACTER*80 TITLE(4)
C

      CHARACTER*10 POTFLG

      CHARACTER*10 LABEL

      REAL*8 EC
      

      REAL*8 VCPRM
C      INTEGER IDPMA
C      INTEGER IDPTR
      CHARACTER*80 FILNAM
      INTEGER I1,I2

      REAL*8 DIPACP(53),
     1       DIPACQ(53),
     2       DIPBCX(53)

c ---- Quantum number parameters for stretching wave functions ----- 
      include 'dipolsys.h'

      COMMON /DIPAC/ DIPACP,DIPACQ
      COMMON /DIPBC/ DIPBCX 

      INTEGER ZC1,ZC2,ZC3,ZELEC,I0,NSTINT
      COMMON /CHARGE/ EC,ZC1,ZC2,ZC3,ZELEC


C      REAL*8 COEFK1,COEFK3

C      COMMON /AMPDI/ IDPMA,IDPTR

C

c      include 'dimen.h'
      include 'molcul.h'
      include 'mu.h'
      include 'rensys.h'



C---------------------------------------------------C
C---------------------------------------------------C
C ------------ START OF THE SUBROUTINE -------------C

C
C     Belegen der Elementarladung EC (kopiert von g.o.)
C
C
C ***** ELEMENTARY CHARGE IN COULOMBS
C
      EC=1.60217733D-19
C
C ***** VACUUM PERMITTIVITY IN FARAD/METER
C
      VCPRM=8.854187817D-12
C
C ***** ELEMENTARY CHARGE IN ELECTROSTATIC UNITS (E.S.U.)
C
      EC=EC/DSQRT(4.0D-9*PI*VCPRM)
C
C ***** ELEMENTARY CHARGE IN DEBYE/ANGSTROM =
C
C       1.0E8 * DEBYE/CM = 1.0E-10 * E.S.U. * CM/CM
C
      EC=1.0D10*EC
C************************************************************

C     ------- OPEN FILES --------- 

      InpUnit = 5
      OutUnit = 7
      TMPUnit = 9
c      FILNAM=DirInput//FileInput
c      OPEN(UNIT=InpUnit,FILE=FILNAM)

c      FILNAM=DirInput//FileOutput
c      OPEN(UNIT=OutUnit,FILE=FILNAM)


C---  SKIP 1 LINE 
      READ (InpUnit,5000) LABEL
      DO I0=1,4
         READ  (InpUnit,5000) TITLE(I0)
      ENDDO
C---  SKIP 7 LINES 
      DO I0=1,7
         READ  (InpUnit,5000) LABEL 
      ENDDO

C---  READ DATA --- 
      READ(InpUnit,5050) DirInput
      READ(InpUnit,5050) FileOutput

      I0=1
      DO WHILE (DirInput(I0:I0).NE.' ')
         I0=I0+1
      ENDDO
      NumDirLetters = I0-1

      I0=1
      DO WHILE (FileOutput(I0:I0).NE.' ')
         I0=I0+1
      ENDDO
      NumFileLetters = I0-1


      FILNAM=DirInput(1:NumDirLetters)//FileOutput(1:NumFileLetters)//
     c  ".out"
      OPEN(UNIT=OutUnit,FILE=FILNAM,status='REPLACE')

      FILNAM=DirInput(1:NumDirLetters)//FileTMP
      OPEN(UNIT=TMPUnit,FILE=FILNAM)

      WRITE (OutUnit,4990)

C     --- WHAT TIME IS IT NOW?
      StartTime = 0
      CALL CLOCKV
      WRITE (OutUnit,5100) (TITLE(I0), I0=1,4)
     
C---  SKIP 5 LINES 
      DO I0=1,5
         READ  (InpUnit,5000) LABEL 
      ENDDO

c------ Calculate all over matrix elements again or use calculated before   ---c
c------ 1 - calculate, 2 - read and use 
        READ (InpUnit,*) CalcOrRead

C---  SKIP 4 LINES 
      DO I0=1,4
         READ  (InpUnit,5000) LABEL 
      ENDDO

C---  READ DATA --- 
C       Jmax  V2max  Na  Nb   (V1+V3)mx 
        READ (InpUnit,*) 
     1        KMAXAB,NV2Max(1),Nvib(1,1),Nvib(1,2),NV1V3Max(1)
        READ (InpUnit,*) 
     1        KMAXAB,NV2Max(2),Nvib(2,1),Nvib(2,2),NV1V3Max(2)
        READ (InpUnit,*) 
     1        KMAXC,NV2MAXC,   Nvib(3,1),Nvib(3,2),NV1V3Max(3)

c------ FINDING MAXIMUM VALUE OF NV1V3Max(i)
        NV2MaxAB=MAX(NV2Max(1),NV2Max(2))

c------ FINDING MAXIMUM VALUE OF NV1V3Max(i)
        MaxV1V3=MAX(NV1V3Max(1),NV1V3Max(2),NV1V3Max(3))

c------ FINDING MAXIMUM VALUE OF NV1V3Max(i)
        MaxNvibMax=MAX(Nvib(1,1),Nvib(1,2),Nvib(2,1),Nvib(2,2),
     1                 Nvib(3,1),Nvib(3,2))



c------ CALCULATING NUMBER OF ALL DIFFER. COMBINATION  (V1,MaxNumberV1V3-V1)
        i0 = 0
        DO I1 = 0, MaxV1V3
        DO I2 = 0, I1
           i0 = i0+1
           ENDDO
        ENDDO
        MaxNumberV1V3 = i0

C---  SKIP 4 LINES
      DO I0=1,4
         READ (InpUnit,5000) LABEL 
      ENDDO

c------ Interval of bands under consideration ---c
        READ (InpUnit,*) MinJC,MaxJC
        READ (InpUnit,*) MinNVibC,MaxNVibC
        READ (InpUnit,*) MinGamVibC,MaxGamVibC
        READ (InpUnit,*) MinV2C,MaxV2C
        READ (InpUnit,*) MinKC,MaxKC
        READ (InpUnit,*) EnerMaxC

C---  SKIP 1 LINES
        READ (InpUnit,5000) LABEL 

        READ (InpUnit,*) MinSurf,MaxSurf
        READ (InpUnit,*) MinJAB,MaxJAB
        READ (InpUnit,*) MinNVibAB,MaxNVibAB
        READ (InpUnit,*) MinGamVibAB,MaxGamVibAB
        READ (InpUnit,*) MinV2AB,MaxV2AB
        READ (InpUnit,*) MinKAB,MaxKAB
        READ (InpUnit,*) EnerMaxAB


C---  SKIP 3 LINES
      DO I0=1,3
         READ (InpUnit,5000) LABEL 
      ENDDO

C---  READ DATA --- 
      READ (InpUnit,*) I0

      IF ( I0.EQ.1 ) THEN
              IntAct = .TRUE.
      ELSEIF( I0.EQ.2 ) THEN 
              IntAct = .FALSE.
      ELSE
              WRITE (OutUnit,9003) I0
              STOP
      ENDIF


      READ (InpUnit,*) AbsTemper
      READ (InpUnit,*) NuMin,NuMax
      READ (InpUnit,*) StrengthLimit
      READ (InpUnit,*) IntensLimit

C---  SKIP 3 LINES
      DO I0=1,3
         READ (InpUnit,5000) LABEL 
      ENDDO

C---  READ DATA --- 
        READ (InpUnit,*) RHOMAX
        READ (InpUnit,*) NSTINT
        RhoMaximum = RHOMAX
        NumPtsRho  = NSTINT
C--------Set bond length requirements
        READ (InpUnit,*) NumPtsBond
        READ (InpUnit,*) RMin(1)
        READ (InpUnit,*) RMin(2)
        READ (InpUnit,*) RMax(1)
        READ (InpUnit,*) RMax(2)

C---  SKIP 3 LINES
      DO I0=1,3
         READ (InpUnit,5000) LABEL 
      ENDDO

C---  READ DATA --- 
      READ (InpUnit,5720) POTFLG
      IF (POTFLG(1:4) .EQ. 'SYMM') THEN
              SYMM=.TRUE.
      ELSEIF (POTFLG(1:4) .EQ. 'UNSY') THEN
              SYMM=.FALSE.
      ELSE
              WRITE (OutUnit,9002) POTFLG
              STOP
      ENDIF
      READ (InpUnit,*) XMULTI
      READ (InpUnit,*) LambdaAB 
      READ (InpUnit,*) LambdaC

C---  SKIP 1 LINE
      READ (InpUnit,5000) LABEL 

C---  READ DATA --- 
      READ (InpUnit,*) REFM1,REFM2,REFM3
      READ (InpUnit,*) M1,M2,M3 

      READ (InpUnit,*) ZC1,ZC2,ZC3

C---  SKIP 1 LINE
      READ (InpUnit,5000) LABEL 

      READ (InpUnit,*) RE12,RE32
      READ (InpUnit,*) RE0(1,1),RE0(2,1)
      READ (InpUnit,*) RE0(1,2),RE0(2,2)
      READ (InpUnit,*) RE0(1,3),RE0(2,3)
C---  SKIP 1 LINE
      READ (InpUnit,5000) LABEL 

      READ (InpUnit,*) AA1,AA3
      READ (InpUnit,*) ZELEC
      READ (InpUnit,*) Te_C
      READ (InpUnit,*) Te_A

      gns = 0
      !      
      IF (SYMM) THEN
        READ (InpUnit,*) gns(1:4)
        GammaLevel = (/'A1','A2','B2','B1'/) 
      ELSE
        READ (InpUnit,*) gns(1:2)
        GammaLevel = (/'A''','A"','XX','XX'/) 
      ENDIF


C---  SKIP 3 LINES
      DO I0=1,3
         READ (InpUnit,5000) LABEL 
      ENDDO

      READ (InpUnit,*) OUT_COEFF

C---  SKIP 4 LINES
      DO I0=1,4
         READ (InpUnit,5000) LABEL 
      ENDDO
     
      IF (SYMM) THEN
        CALL RDPSYY(InpUnit,DIPACQ)
C---    SKIP 1 LINE
        READ (InpUnit,5000) LABEL 
        CALL RDPSYZ(InpUnit,DIPACP)
C---    SKIP 1 LINE
        READ (InpUnit,5000) LABEL 
        CALL RDPSYY(InpUnit,DIPBCX)
      ELSE  
        CALL RDPAS(InpUnit,DIPACQ)
C---    SKIP 1 LINE
        READ (InpUnit,5000) LABEL 
        CALL RDPAS(InpUnit,DIPACP)
C---    SKIP 1 LINE
        READ (InpUnit,5000) LABEL 
        CALL RDPAS(InpUnit,DIPBCX)
      ENDIF
C--------------------------------------------------C
C--------------------------------------------------C


        IF ( CalcOrRead.EQ.1 ) THEN 
             WRITE(OutUnit,6331)
        ELSE
             WRITE(OutUnit,6332)
        ENDIF 

        WRITE(OutUnit,6329)         
     1          MinSurf,MaxSurf,
     1          MinJAB,MaxJAB,
     3          MinNVibAB,MaxNVibAB,
     4          GammaLevel(2*MinGamVibAB-1),
     5          GammaLevel(2*MaxGamVibAB-1),
     6          MinV2AB,MaxV2AB,
     7          MinKAB,MaxKAB,
     2          MinJC,MaxJC,
     8          MinNVibC,MaxNVibC, 
     9          GammaLevel(2*MinGamVibC-1),
     1          GammaLevel(2*MaxGamVibC-1),
     2          MinV2C,MaxV2C,
     3          MinKC,MaxKC

        WRITE(OutUnit,6330) EnerMaxAB,EnerMaxC,NuMin,Numax

        WRITE(OutUnit,6319) RHOMAX,NumPtsRho,
     1                      KMAXAB,NV2MaxAB,KMAXC,NV2MaxC
        WRITE(OutUnit,6320) NV1V3Max(1),Nvib(1,1),Nvib(1,2),
     1                      Nvib(2,1),Nvib(2,2)
        WRITE(OutUnit,6321) NumPtsBond,RMin(1),RMin(2),RMax(1),RMax(2)

C        MaxNvib = MAX(Nvib(1,1),Nvib(1,2),
C     1                Nvib(2,1),Nvib(2,2))
C

        WRITE (OutUnit,7712) XMULTI
        WRITE (OutUnit,7710) LambdaAB,LambdaC 
  
C
      WRITE (OutUnit,5800) RE12,RE32,
     1                     RE0(1,1),RE0(2,1),
     1                     RE0(1,2),RE0(2,2),
     1                     RE0(1,3),RE0(2,3)


      WRITE (OutUnit,5900) AA1,AA3
C
      WRITE(OutUnit,5709) StrengthLimit,IntensLimit,AbsTemper
      WRITE(OutUnit,5706) ZC1,ZC2,ZC3
      WRITE(OutUnit,5707) 'ZELEC=    ',ZELEC

C
C  Write out dipole parameter values
C
        WRITE(OutUnit,5711)
        CALL WDIP(OutUnit,'Dipole Moment ACQ',DIPACQ)
        CALL WDIP(OutUnit,'Dipole Moment ACP',DIPACP)
        CALL WDIP(OutUnit,'Dipole Moment BCX',DIPBCX)
C     

C

C------------- close input file -------------------
        CLOSE(InpUnit)
C------------- close input file -------------------




      RETURN

C ------------- END OF THE SUBROUTINE --------------C
C---------------------------------------------------C
C---------------------------------------------------C


5305  FORMAT('1     ',A8)
4990  FORMAT('1','ABCTOC.INP.INF  INPUT PROCESSING ROUTINE',
     1      ' CALLED'/)
5000  FORMAT(A80)
5050  FORMAT(A80)
5100  FORMAT('0',3(132('*')/' '),16('*'),100X,16('*')/
     1      4(' ',16('*'),10X,A80,10X,16('*')/),
     2        ' ',16('*'),100X,16('*')/' ',3(132('*')/' '))
5701  FORMAT(A10,3F20.12)
5720  FORMAT(A6)
9002  FORMAT('0','RENTEL.INP.ERR WRONG POTENTIAL TYPE ',A6,' SPECIFIED'/
     1      ' ','                ( SYMM AND UNSYMM ARE ALLOWED)'//)
9003  FORMAT('0','RENTEL.INP.ERR WRONG TYPE OF SPECTRUM ',I2,
     1           ' SPECIFIED'/
     1      ' ','                ( 1 OR 2 )'//)
9004  FORMAT('0','RENTEL.INP.ERR WRONG SPECTRUM DIRECTION ',I2,
     1           ' SPECIFIED'/
     1      ' ','                ( 1 OR 2 )'//)
5800  FORMAT('0',17X,'    REF12 =',F10.6,' ANGSTROM,   REF32 =',
     1      F10.6,' ANGSTROM'/
     1       '0',17X,'A   RE012 =',F10.6,' ANGSTROM,   RE032 =',
     1      F10.6,' ANGSTROM'/
     1       '0',17X,'B   RE012 =',F10.6,' ANGSTROM,   RE032 =',
     1      F10.6,' ANGSTROM'/
     1       '0',17X,'C  RE012 =',F10.6,' ANGSTROM,   RE032 =',
     1      F10.6,' ANGSTROM'/)
5900  FORMAT('0',17X,'    A1 =',F8.4,' ANGSTROM-1,   A3 =',F8.4,
     1      ' ANGSTROM-1'/)
5990  FORMAT('0',5X,12('*'),' POTENTIAL ENERGY EXPANSION COEFFI',
     1      'CIENTS (CM-1) ',12('*')/
     1       '0',17X,'   F11     =',F16.5,'   F33     =',F16.5)
5705  FORMAT(A10,3I4)
5702  FORMAT(A10,I4)
5504  FORMAT('0',5X,12('*'),' Nuclear Spin Statistical ',
     1       'Weights ',12('*'))
5503  FORMAT(A10,2I4)
5703  FORMAT(A10,D16.8)
5704  FORMAT(A10,4I4)
5711  FORMAT('0',5X,12('*'),' Dipole Moment Parameters ',12('*')//)
5706  FORMAT('0',5X,'ZC1= ',I4,' ZC2= ',I4,' ZC3= ',I4)
5707  FORMAT('0',5X,A10,I4)
5709  FORMAT('0',5X,'Output limit for linestrength is ',F12.6/
     1       '0',5X,'Output limit for Intensity    is ',D12.4/
     1       '0',5X,'AbsTemper  = ',F12.5,' (ABSOLUTE TEMPERATURE)')

5710  FORMAT(A10,I5)

6319  FORMAT('0',5X,12('*'),' NUMERICAL PARAMETERS ',
     1  '(BENDING) ',12('*')//,
     2  '0',17X,'  RHOMAX  = ',F10.5,' RADIANS'/
     3  '0',17X,'  NumPtsRho  = ',I5,' (NUMBER OF INTEGRATION STEPS)'/
     4  '0',17X,'  KMAX AB = ',I5,' (MAX NUMBER OF | K  >)'/
     5  '0',17X,'  NV2MAX AB  = ',I5,' (MAX NUMBER OF | V2 >)'/
     4  '0',17X,'  KMAX C  = ',I5,' (MAX NUMBER OF | K  >)'/
     5  '0',17X,'  NV2MAX C  = ',I5,' (MAX NUMBER OF | V2 >)'/)

6327  FORMAT('0',5X,12('*'),' EQUALIBRIUM ENERGY '/
     1  '0',17X,'  Te A state = ',F10.2,' 1/CM '/
     2  '0',17X,'  Te C state = ',F10.2,' 1/CM ')
6328  FORMAT('0',5X,12('*'),' ZERO POINT ENERGY '/
     1  '0',17X,'  ZeroEnergy A state = ',F10.2,' 1/CM '/
     2  '0',17X,'  ZeroEnergy C state = ',F10.2,' 1/CM ')


6320  FORMAT('0',5X,12('*'),' NUMERICAL PARAMETERS ',
     1  '(STRETCHING) ',12('*')//,
     2  '0',17X,'A V1V3Max  = ',I5,' (Max number of |V1+V3>)'/
     3  '0',17X,'B V1V3Max  = ',I5,' (Max number of |V1+V3>)'/
     4  '0',17X,'C V1V3Max  = ',I5,' (Max number of |V1+V3>)'/
     5  '0',17X,'A NumberA1 = ',I5,
     6                    ' (Number of stretch. fncts A1 symm)'/
     7  '0',17X,'A NumberB2 = ',I5,
     8                    ' (Number of stretch. fncts B2 symm)'/
     9  '0',17X,'B NumberA1 = ',I5,
     1                    ' (Number of stretch. fncts A1 symm)'/
     2  '0',17X,'B NumberB2 = ',I5,
     3                    ' (Number of stretch. fncts B2 symm)'/
     4  '0',17X,'C NumberA1 = ',I5,
     5                    ' (Number of stretch. fncts A1 symm)'/
     6  '0',17X,'C NumberB2 = ',I5,
     7                    ' (Number of stretch. fncts B2 symm)'/
     8  '0',17X,'  NvibA = ',I5,
     9                    ' (Max number of N vib for A state)'/
     1  '0',17X,'  NvibB = ',I5,' (                   for B state)'/
     2  '0',17X,'  NvibC = ',I5,' (                  for C state)')
6321  FORMAT(
     1  '0',17X,'NumPtsBond = ',I5,' (NUMBER OF INTEGRATION STEPS)'/
     2  '0',17X,'R12Min     = ',F10.5,' (Min value of R12 bond length)'/
     2  '0',17X,'R32Min     = ',F10.5,' (Min value of R32 bond length)'/
     2  '0',17X,'R12Max     = ',F10.5,' (Max value of R12 bond length)'/
     2  '0',17X,'R32Max     = ',F10.5,' (Max value of R32 bond length)')
6329  FORMAT('0',5X,12('*'),' CONSIDER TRANSITIONS ONLY BETWEEN BANDS'/
     1  '0',17X,'FROM AB STATE eta = [',I2,'..',I2,'] ',
     1  '0',17X,'FROM AB STATE J = [',I2,'..',I2,'] ',
     1          'AND NVib = [',I2,'..',I2,'] ',
     1          'AND GamVib = [',A2,'..',A2,'] ', 
     1          'AND V2 = [',I2,'..',I2,']', 
     1          'AND Ka = [',I2,'..',I2,']'/
     1  '0',17X,'TO   C  STATE J = [',I2,'..',I2,'] ',
     1          'AND  NVib = [',I2,'..',I2,'] ',
     1          'AND GamVib = [',A2,'..',A2,'] ', 
     1          'AND V2 = [',I2,'..',I2,']', 
     1          'AND Ka = [',I2,'..',I2,']')
6330  FORMAT('0',5X,12('*'),' CONSIDER ALL TRANSITIONS FROM A-B STATE',
     1                      ' LOWER THEN',F10.2,' 1/CM'/
     1       '0',5X,12(' '),'                          TO   C   STATE',
     1                      ' LOWER THEN',F10.2,' 1/CM'/
     1       '0',5X,12(' '),' IN INTERVAL [',F10.2,'...',F10.2,']') 



6331  FORMAT('0',5X,12('*'),' IN THIS SATION WE CALCULATE ALL OVER'
     1                      ' MATRIX ELEMENTS AGAIN')
6332  FORMAT('0',5X,12('*'),' IN THIS SATION WE READ AND USE'
     1                      ' CALCULATED MATRIX ELEMENTS')

7712  FORMAT('0',5X,12('*'),' SPIN-ORBIT COUPLING PARAMETERS',
     1      12('*')//
     2      17X,'   XMULTI   =',F16.5)



7710  FORMAT('0',5X,12('*'),' ELECTRONIC ANGULAR MOMENTUM PAR',
     1      'AMETER  ',12('*')//
     2      17X,'   XLAMB    =',F16.5)



      END

C###################      ReadInputData        #######################C
C#####################################################################C



C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : 
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : CHECK INPUT DATA 
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################
      SUBROUTINE CheckInputData 
C     IMPLICIT REAL*8 (A-H,O-Z)
C     IMPLICIT LOGICAL (A-Z)

      include 'mu.h'
      include 'dipolsys.h'
      include 'molcul.h'

      LOGICAL File_exist
      CHARACTER*80 FileName(40),FilNam
      INTEGER k0,MaxBounds,Num1,Num2,Num3
      REAL*8 Val1,Val2,NearZero
      CHARACTER*2 NumChar

C ----------------------------------------------------------------------C
C ----------------------------------------------------------------------C


      NearZero = 1.0D-6 



C     --- CHECKING OF THE EXISTING OF FILES 


      
      FileName(1)=DirInput(1:NumDirLetters)//FileStretchAB
      FileName(2)=DirInput(1:NumDirLetters)//FileStretchC
      FileName(3)=DirInput(1:NumDirLetters)//FileGamma   
      FileName(4)=DirInput(1:NumDirLetters)//FileBendAB
      FileName(5)=DirInput(1:NumDirLetters)//FileBendC   
      FileName(6)=DirInput(1:NumDirLetters)//FileCoefAB
      FileName(7)=DirInput(1:NumDirLetters)//FileCoefC



      DO K0 = 1,7
        FilNam = FileName(K0)
        INQUIRE(FILE=FilNam,EXIST=File_exist)
        IF ( .NOT.File_exist ) GOTO 91
      ENDDO 

      
      DO K0 = 1,2
        FilNam = FileName(K0)
        OPEN(UNIT=92,FILE=FilNam)
        READ(92,*) Val1,Val2
        IF ( DABS(Val1-AA1).GT.1E-4 .OR.
     1       DABS(Val2-AA3).GT.1E-4 ) GOTO 92
        READ(92,*) Val1,Val2
C        IF ( DABS(Val1-Kappa(1,2*K0-1)) .GT.NearZero .OR.
C     1       DABS(Val2-Kappa(2,2*K0-1)) .GT.NearZero  ) GOTO 92
        READ(92,*) Val1,Val2
c        IF ( DABS(Val1-RE12).GT.0.001 .OR. 
c     1       DABS(Val2-RE32).GT.0.001 ) GOTO 92
        READ(92,*) Num1
        READ(92,*) Num1,Num2
        IF ( Num1.NE.Nvib(2*K0-1,1) .OR. 
     1       Num2.NE.Nvib(2*K0-1,2) ) GOTO 92
      ENDDO 


C      IF (   DABS(Kappa(1,1)-Kappa(1,2)).GT.NearZero .OR.
C     1       DABS(Kappa(2,1)-Kappa(2,2)).GT.NearZero ) GOTO 92


        FilNam = FileName(3)
        OPEN(UNIT=92,FILE=FilNam)

        READ(92,*) Num1,Num2
        IF ( Num1.NE.NumPtsRho .OR. Num2.NE.KMAXAB+1) GOTO 92
        READ(92,*) Val1,Val2
        IF ( DABS(Val1-RhoMaximum/NumPtsRho).GT.NearZero .OR. 
     1       DABS(Val2-RhoMaximum).GT.NearZero ) GOTO 92

        FilNam = FileName(4)
        OPEN(UNIT=92,FILE=FilNam)
        READ(92,*) Num1,Num2,Num3

        IF ( Num1.NE.NumPtsRho .OR. Num2.NE.KMAXAB+1 .OR. 
     1       Num3.NE.NV2Max(1)+1) GOTO 92

        READ(92,*) Val1,Val2
        IF ( ABS(Val1-RhoMaximum/NumPtsRho).GT.NearZero .OR. 
     1       ABS(Val2-RhoMaximum).GT.NearZero ) GOTO 92

        FilNam = FileName(5)
        OPEN(UNIT=92,FILE=FilNam)
        READ(92,*) Num1,Num2,Num3
        IF ( Num1.NE.NumPtsRho .OR. Num2.NE.KMAXC+1 .OR. 
     1       Num3.NE.NV2MaxC+1) GOTO 92
        READ(92,*) Val1,Val2
        IF ( DABS(Val1-RhoMaximum/NumPtsRho).GT.NearZero .OR. 
     1       DABS(Val2-RhoMaximum).GT.NearZero ) GOTO 92


      GOTO 93 
      
c----------------------------c
c---------  STOP ------------c
92       CLOSE(92) 
              WRITE(*,*) 'ERROR: INPUT PARAMETERS FROM FILE', 
     1   FilNam,'DO NOT AGREE WITH PARAMETERS FROM THE INPUT FILE'
         WRITE (*,*) VAL1,VAL2,NUM1,NUM2
c          PAUSE '-- PRESS ANY KEY --'
          STOP
c---------  STOP ------------c
c----------------------------c

93       CLOSE(92) 

C     --- CHECKING OF THE BOUNDS OF ARRAYS 


      IF (  MaxNvibMax.GT.UpperNvib ) THEN 
         WRITE(*,*) 'ERROR: PARAMETER UpperNvib', 
     1  ' TOO SMALL TO BE A BOUND OF THE ARRAY ',
     1  ' NotateStretch and others, LESS THAN  MaxNvibMax = ',MaxNvibMax
          PAUSE '-- PRESS ANY KEY --'
          STOP
      ENDIF


      IF (  NV2MaxAB.GT.UpperV2 .OR. NV2MaxC.GT.UpperV2)  THEN 
         WRITE(*,*) 'ERROR: PARAMETER UpperV2', 
     1  ' TOO SMALL TO BE A BOUND OF THE ARRAY ',
     1  ' LESS THAN  NV2Max = ',NV2MaxAB,' OR ',NV2MaxC
          PAUSE '-- PRESS ANY KEY --'
          STOP
      ENDIF

      IF (  KMAXAB.GT.UpperK .OR. KMAXC.GT.UpperK)  THEN 
         WRITE(*,*) 'ERROR: PARAMETER UpperK', 
     1  ' TOO SMALL TO BE A BOUND OF THE ARRAY ',
     1  ' LESS THAN  KMAX = ',KMAXAB,' OR ',KMAXC
          PAUSE '-- PRESS ANY KEY --'
          STOP
      ENDIF


      IF (  NumPtsRho.GT.MaxNumRhoPoints ) THEN 
         WRITE(*,*) 'ERROR: PARAMETER MaxNumRhoPoints', 
     1  ' TOO SMALL TO BE A BOUND OF THE ARRAY ',
     1  ' CoefMuRho, LESS THAN  NumPtsRho = ',NumPtsRho            
          PAUSE '-- PRESS ANY KEY --'
          STOP
      ENDIF


      IF (  NumPtsBond.GT.MaxNumRhoPoints ) THEN 
         WRITE(*,*) 'ERROR: PARAMETER MaxNumRhoPoints', 
     1  ' TOO SMALL TO BE A BOUND OF THE ARRAY ',
     1  ' FnctsMorse, LESS THAN  NumPtsBond = ',NumPtsBond
          PAUSE '-- PRESS ANY KEY --'
          STOP
      ENDIF


c      IF (  MOD(NumPtsBond,2).EQ.1 ) THEN 
c               WRITE(*,*) 'ERROR: PARAMETER NumPtsBond', 
c     1  ' MUST BE EVEN! '
c                  PAUSE '-- PRESS ANY KEY --'
c                  STOP
c      ENDIF
c
c
c      IF (  MOD(NumPtsRho,2).EQ.1 ) THEN 
c         WRITE(*,*) 'ERROR: PARAMETER NumPtsRho', 
c     1  ' MUST BE EVEN! '
c          PAUSE '-- PRESS ANY KEY --'
c          STOP
c      ENDIF




C     ------------- END OF CHECKING -------------- C
      GOTO 90    !  --- NO PROBLEM! 


91    WRITE(*,*) 'ERROR: INPUT FILE '//FilNam//' DOES EXIST!'
      PAUSE '-- PRESS ANY KEY --'
      STOP


90    CONTINUE 

      RETURN

      END

C###################       CheckInputData      #######################C
C#####################################################################C




C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : 
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : CLOCK
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################

      SUBROUTINE CLOCKV
C      REAL*8 NEWTIM
C      REAL XI1,XARRAY(2)

      CHARACTER*8 HOUR
      CHARACTER*2 CH_SEC,CH_MIN,CH_HOUR
      INTEGER NTemp
      INTEGER Hour1,Min1,Sec1,TimePassedToday,Date1
      INTEGER SystemTime,CR,CM
      include 'mu.h'

      CALL SYSTEM_CLOCK(SystemTime,CR,CM)

      TimePassed = INT(SystemTime/CR)

      IF ( SystemTime.EQ.0 ) THEN 
          CALL TIME (HOUR)
          CH_SEC= HOUR(7:)
          CH_MIN= HOUR(4:5)
          CH_HOUR=HOUR(1:2)

          READ (CH_SEC ,*)  SEC1
          READ (CH_MIN ,*)  MIN1
          READ (CH_HOUR,*) HOUR1

          WRITE(6,1124) Hour1,Min1,Sec1
      ELSE 

      Date1 = INT(TimePassed/86400)
      NTemp = TimePassed - Date1*86400
      Hour1 = INT(NTemp/3600)
      NTemp = NTemp - Hour1*3600
      Min1  = INT(NTemp/60)
      Sec1 = NTemp - Min1*60

C      WRITE(OutUnit,1125) Hour1,Min1,Sec1

      ENDIF 


1124  FORMAT(1H0,' ABTOC.INT.INF   START TIME: ',I3,':',I3,':',I3 /)
1125  FORMAT(1H0,' ABTOC.INT.INF   ',I3,':',I3,':',I3)
3001  FORMAT(1H0,' ABTOC.INT.INF   START INTENSITY CALCULATIONS '/)
      RETURN
      END
C###################          CLOCKV           #######################C
C#####################################################################C

