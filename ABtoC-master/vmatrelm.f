!!!*************************************************************
! 文件/File: vmatrelm.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vmatrelm.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


C##########################################################################
C#             :
C# DATE        : 05.06.2000
C# AUTHOR      : Yurchenko
C# UPDATES     :
C# LANGUAGE    : FORTRAN
C# PURPOSE     : CALCULATION OF THE DIPOLE MOMENT MATRIX ELEMENTS ON THE 
C#             : SYMMETRIZED STRETCHING FUNCTIONS |NVibGammaVib> 
C#             : AND BENDING FUNCTIONS |V2 K> 
C#
C###########################################################################

      SUBROUTINE  Calc_VibrationMatrixElement

C      IMPLICIT REAL*8 (A-H,O-Z)
C     IMPLICIT LOGICAL (A-Z)
C
      include 'dipolsys.h'
      include 'mu.h'
C
C      REAL*8 VibrationMatrixElement(MaxNumVibMatEl)
C
C      COMMON /VibrMatEl/ VibrationMatrixElement
C
C      INTEGER IndexME(2,MaxNumVibMatEl),
C     1        MaxTermK(MaxNumJ+1,2),MinTermK(2)
C
C      COMMON /IndexMatEl/ IndexME,MaxTermK,MinTermK

      REAL*8  StretchMatrixElements(15,2,2,2,UpperNvib,UpperNvib)
      COMMON /DIPSTRETCH/ StretchMatrixElements


      REAL*8 DipMomentBend(2,2,UpperV2,UpperV2,UpperK,2,15)
      COMMON /DIPBEND/ DipMomentBend


      
      LOGICAL File_exist
      CHARACTER*80 FILNAM


C--------------------------------------------------------------------C
C--------------------------------------------------------------------C
C --------------------- START OF THE SUBROUTINE ---------------------C


      WRITE(6,3001)
3001  FORMAT(1H0,' ABTOC.INT.INF   START MATRIX ELEMENTS CALCULATIONS'/)


C     --- CHECKING OF THE EXISTING FILES WITH DIPOLE MOMENT 
      FILNAM=DirInput(1:NumDirLetters)//FileVibME
      INQUIRE(FILE=FilNam,EXIST=File_exist)


C     ---- CHECK NECESSITY OF THE CALCULATION 
      IF ( CalcOrRead.EQ.1 .OR. .NOT. File_exist ) THEN 


C     ---- CALCULATE MATRIX ELEMENTS OF THE DIPOLE MOMENT 
C     ---- ON THE STRETCHING WAVE FUNCTIONS 

      CALL StretchMatrixElement 

C     ---- CALCULATE MATRIX ELEMENTS OF THE DIPOLE MOMENT OF THE 
C     ---- TRIATOMIC MOLECULE USING A NUMERICAL INTAGRATION 
      CALL BendMatrixElement


C     ---- CALCULATION OF THE DIPOLE MOMENT MATRIX ELEMENTS 
C     ---- SEPARATLY FOR X,Y OR Z COMPONENT AND FOR COS OR SIN WEIGHT 

      CALL   DipolMomentVibMatrElem


      ENDIF 


      RETURN
      END

C##################       VibrationMatrixElement       #####################
C###########################################################################


C##########################################################################
C#             :
C# DATE        : 05.06.2000
C# AUTHOR      : Yurchenko
C# UPDATES     :
C# LANGUAGE    : FORTRAN
C# PURPOSE     : CALCULATION OF THE DIPOLE MOMENT MATRIX ELEMENTS ON THE 
C#             : SYMMETRIZED STRETCHING FUNCTIONS |NVibGammaVib> 
C#             : AND BENDING FUNCTIONS |V2 K> SEPARATLY FOR X,Y OR Z COMPONENT  
C#             : AND FOR COS OR SIN WEIGHT 

      SUBROUTINE  DipolMomentVibMatrElem

C             : 
C
C    I  Used Symm     
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


C###########################################################################

      include 'dipolsys.h'
      include 'mu.h'
C
C      REAL*8  VibrationMatrixElement(MaxNumVibMatEl)
      INTEGER*4 IndexTerm_AB,IndexTerm_C 

      REAL*8  StretchMatrixElements(15,2,2,2,UpperNvib,UpperNvib)
      COMMON /DIPSTRETCH/ StretchMatrixElements

      REAL*8 DipMomentBend(2,2,UpperV2,UpperV2,UpperK,2,15)


      COMMON /DIPBEND/ DipMomentBend

      INTEGER  sigma,IPowers, ISurfL, GammaVibL, GammaVibR, 
     1         NVibL, NVibR, KL, KR, V2L, V2R, Term,
     1         MinTermK(2)



      INTEGER  MaxTerm,NDeltaK,DeltaK,theta

      REAL*8 MatrElemTemp

      CHARACTER*80 FILNAM


C--------------------------------------------------------------------C
C--------------------------------------------------------------------C
C --------------------- START OF THE SUBROUTINE ---------------------C



      WRITE(6,3001)
3001  FORMAT(1H0,' ABTOC.INT.INF   FINAL MATRIX ELEMENTS '/)


C*********************************************************************
C
C     OPEN FILE FOR STORING THE DIPOLE MOMENT MATRIX ELEMENTS 
C
C*********************************************************************
      FILNAM=DirInput(1:NumDirLetters)//FileVibME
      OPEN(UNIT=31,FILE=FILNAM)
      IF ( CalcOrRead.EQ.2 )  WRITE(6,500) FilNam

      Term = 0

C     LOOP OVER DIFFERENT QUANTUM NUMBERS OF THE SYMMTR. STRETCH. EIGEN. FUNCT. OF THE H_str 
C

      DO theta = 0,1                 ! --- theta
      DO KL =0,KMAXAB                ! -- KLeft   - bend quantum value 
      DO DeltaK = -theta,theta,2     ! --- DeltaK 
           KR = KL + DeltaK

      IF ( KR.GE.0 .AND. KR.LE.KMAXC ) THEN     ! --- KR <=KMAXC 
C     ---- NDeltaK  - index for distinguishing DeltaK = -1,0,1 => NdeltaK = 1,1,2
           NDeltaK = INT((DeltaK+2)/3)+1

      DO ISurfL = 1,2                  ! --- A,B - electronic states  Left  < |
      DO V2L=0,NV2MAX(ISurfL)          ! -- V2Left  - bend quantum value 
      DO V2R=0,NV2MAXC                 ! -- V2Right - bend quantum value 
      DO GammaVibL=1,2                          ! -- A1,B2 Left  < |
      DO GammaVibR=1,2                          ! -- A1,B2 Right | >
      DO NVibL = 1,Nvib(ISurfL,GammaVibL)  ! -- NVibLeft < | 
      DO NVibR = 1,Nvib(3,GammaVibR)  ! -- NVibRight| >



C    
C     LOOP OVER DIFFERENT POWERS OF THE Y_1^i*Y_3^j
C
      
      MatrElemTemp = 0.0D+00
      DO IPowers = 1,15             ! ---- Y Powers
         MatrElemTemp = MatrElemTemp + 
     1   DipMomentBend(theta+1,ISurfL,V2L+1,V2R+1,
     1                   KL+1,NDeltaK,IPowers)*
     2   StretchMatrixElements(IPowers,ISurfL,
     3                         GammaVibL,GammaVibR,NVibL,NVibR)
	   continue 
      ENDDO             ! ---- YPowers


      IF ( DABS(MatrElemTemp).GT.AlmostZero  ) then
     
             Term = Term + 1

c      IF ( GammaVibL.NE.GammaVibR  .AND. lambda.EQ.2 ) THEN   
c              PAUSE '-- GammaVibL  = GammaVibR --'
c     ENDIF 


          IF ( Term .GT. MaxNumVibMatEl ) THEN
            WRITE(*,*) 'ERROR: PARAMETER MaxNumVibMatEl', 
     1      ' TOO SMALL TO BE A BOUND OF THE  ARRAY',
     1      ' VibrationMatrixElement, LESS THAN ',Term
              STOP
          ENDIF

c----  Index matrix for left AB state
             IndexTerm_AB = 10000000*ISurfL+
     1                      100000*KL+1000*V2L+10*NVibL+GammaVibL

c----  Index matrix for right C state
             IndexTerm_C  = 100000*KR+1000*V2R+10*NVibR+GammaVibR


c         IndexTerm_AB = (GammaVibL-1)*2+(ISurfL-1)
c         IndexTerm_C  = (GammaVibR -1)*2
c
c         CALL MVBITS(KL   ,0,4,IndexTerm_AB,2)
c         CALL MVBITS(V2L   ,0,4,IndexTerm_AB,6)
c         CALL MVBITS(NVibL,0,4,IndexTerm_AB,10)
c
c         CALL MVBITS(KR   ,0,4,IndexTerm_C,2)
c         CALL MVBITS(V2R  ,0,4,IndexTerm_C,6)
c         CALL MVBITS(NVibR,0,4,IndexTerm_C,10)


             WRITE (31,505) MatrElemTemp,
     1                      IndexTerm_AB,IndexTerm_C


      ENDIF             ! --  MatrElemTemp >0


      ENDDO             ! --  NVibRight| >
      ENDDO             ! --  NVibLeft < | 
      ENDDO             ! --  A1,B2 Right | >
      ENDDO             ! --  A1,B2 Left  < |
      ENDDO             ! --- A,B - electronic states  Left  < |
      ENDDO             ! --  V2Right - bend quantum value 
      ENDDO             ! --  V2Left  - bend quantum value 
      ENDIF             ! --  KR <=KMAXC 
      ENDDO             ! --  DeltaK = KR-KL 

      ENDDO             ! --  KLeft   - bend quantum value 
      ENDDO             ! --- theta


C     --- CHECKING OF THE BOUNDS OF ARRAYS 


      MaxTerm = Term 

      IF ( MaxTerm .GT. MaxNumVibMatEl ) THEN
        WRITE(*,*) 'ERROR: PARAMETER MaxNumVibMatEl', 
     1  ' TOO SMALL TO BE A BOUND OF THE  ARRAY',
     1  ' VibrationMatrixElement, LESS THAN ',MaxTerm
          PAUSE '-- PRESS ANY KEY --'
          STOP
      ENDIF


      CLOSE(31)

      RETURN


505   FORMAT (D24.12,I10,I10)
500   FORMAT (1H0,' ABTOC.INT.INF  : INPUT FILE ',A30,' DOES EXIST!'/
     1'                     : MATRIX ELEMENTS HAVE NOT BEEN CALCULATED'/
     1        '                     : AND ARE COMPUTING NOW')

      END

C#################     DipolMomentVibMatrElem     ####################
C###########################################################################







C##########################################################################
C#             :
C# DATE        : 05.06.2000
C# AUTHOR      : Yurchenko
C# UPDATES     :
C# LANGUAGE    : FORTRAN
C# PURPOSE     : READ OF THE DIPOLE MOMENT MATRIX ELEMENTS ON THE 
C#             : SYMMETRIZED STRETCHING FUNCTIONS |NVibGammaVib> 
C#             : AND BENDING FUNCTIONS |V2 K> SEPARATLY FOR X,Y OR Z COMPONENT  
C#             : AND FOR COS OR SIN WEIGHT 
C###########################################################################


      SUBROUTINE  Read_VibrationMatrixElement

C
      include 'dipolsys.h'
      include 'mu.h'
C
      REAL*8 VibrationMatrixElement(MaxNumVibMatEl)

      COMMON /VibrMatEl/ VibrationMatrixElement

      INTEGER IndexME(2,MaxNumVibMatEl),
     1        MaxTermK(MaxNumJ+1,2),MinTermK(2)

      COMMON /IndexMatEl/ IndexME,MaxTermK,MinTermK

      INTEGER  KL, KR, Term

      INTEGER  MaxTerm,theta,KLOLD,etaAB,Ntemp

      CHARACTER*80 FILNAM


C--------------------------------------------------------------------C
C--------------------------------------------------------------------C
C --------------------- START OF THE SUBROUTINE ---------------------C


      WRITE(6,3001)
3001  FORMAT(1H0,' ABTOC.INT.INF   READ MATRIX ELEMENTS FROM FILE'/ )


C*********************************************************************
C
C     OPEN FILE FOR STORING THE DIPOLE MOMENT MATRIX ELEMENTS 
C
C*********************************************************************
      FILNAM=DirInput(1:NumDirLetters)//FileVibME
      OPEN(UNIT=31,FILE=FILNAM,STATUS='OLD')

      Term = 0

      theta = 0 
      MinTermK(theta+1) = Term +1 
      KLOLD = 0
      MaxTermK(KLold+1,theta+1) = 0

      DO WHILE (.NOT. EOF(31))

         Term = Term + 1
    
         READ(31,505) VibrationMatrixElement(Term),
     1                IndexME(1,Term),
     1                IndexME(2,Term)

c         CALL MVBITS(IndexME(1,Term), 2,4,KL,0)
c         CALL MVBITS(IndexME(2,Term), 2,4,KR,0)

         etaAB=  INT(IndexME(1,Term)/1000000)
         Ntemp = IndexME(1,Term)-etaAB*1000000
         KL=   INT(Ntemp/100000)
         KR=   INT(IndexME(2,Term)/100000)


C         KL=INT(IndexME(1,Term)/10000000)
C         KR=INT(IndexME(2,Term)/10000000)

         IF ( KL.NE.KR .AND. theta.EQ.0 ) THEN 
            MaxTermK(KLold+1,theta+1) = Term - 1
            KLOLD = KL 
            theta = 1 
            MinTermK(theta+1) = Term
         ELSEIF ( KL.NE.KLOLD ) THEN 
            MaxTermK(KLold+1,theta+1) = Term - 1
            KLOLD = KL 
         ENDIF 


C     --- CHECKING OF THE BOUNDS OF ARRAYS 
      IF ( term .GT. MaxNumVibMatEl ) THEN
        WRITE(*,*) 'ERROR: PARAMETER MaxNumVibMatEl', 
     1  ' TOO SMALL TO BE A BOUND OF THE  ARRAY',
     1  ' VibrationMatrixElement, LESS THAN ',Term
          STOP
      ENDIF

      ENDDO            ! --- term 

      MaxTermK(KL+1,theta+1) = Term


      CLOSE(31)

      RETURN


505   FORMAT (D24.12,I10,I10)
      END

C##############         Read_VibrationMatrixElement       ##################
C###########################################################################

