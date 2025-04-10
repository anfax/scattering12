!!!*************************************************************
! 文件/File: intens.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: intens.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      THIS FILE INCLUDES SUBROUTINES DEVOTED TO CALCULATION OF THE C
C      INTENSITIES OF THE ABA THREEATOMIC MOLECULES                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : CALCULATES LINESTRENGTH OF THE TRIATOMIC MOLECULE
C SUBPROGRAMS :
C##########################################################################
      SUBROUTINE Calc_LineStrength

C     IMPLICIT REAL*8 (A-H,O-Z)
C     IMPLICIT LOGICAL (A-Z)

      include 'dipolsys.h'
      include 'mu.h'
      include 'rensys.h'
      include 'molcul.h'

      REAL*8 VibrationMatrixElement(MaxNumVibMatEl)

      COMMON /VibrMatEl/ VibrationMatrixElement

      INTEGER IndexME(2,MaxNumVibMatEl),
     1        MaxTermK(MaxNumJ+1,2),MinTermK(2)

      REAL*8 ThreeJsymbol(KMAXAB+1,3,KMAXAB+1,3),
     1       SixJsymbol(KMAXAB+1,KMAXAB+1,3,3)

      COMMON /IndexMatEl/ IndexME,MaxTermK,MinTermK

      REAL*8 CoeffAB(0:MaxNumCoeffAB),CoeffC(0:MaxNumCoeffC)
      INTEGER IndexAB(2,KMAXAB+1,KMAXAB+1,NV2MAXAB+1,MaxNvibMax,2),
     1        IndexC            (KMAXC+1,NV2MAXC +1,MaxNvibMax,2)

      REAL*8 STemp,RotMatrElem,
     1       MatrSigma,Q_AB,Q_C,AbsInt,Q

      REAL*8 EnergyC(MaxNumLevelC),EnergyAB(MaxNumLevelAB)

      CHARACTER*80 FILNAMAB,FILNAMC

      INTEGER K_AB,Ntemp,V2AB,NVibAB,GammaVibAB,etaAB,
     1        TermAB,TermC,K_C,V2C,NVibC,GammaVibC,
     1        Term,
     1        NotateLevelAB(MaxNumLevelAB,11),
     2        NotateLevelC(MaxNumLevelC,11)

      INTEGER J2_f,K_f,V2_f,V1_f,V3_f,NVib_f,Gamma_f,
     1        tau_f,UNIT_f,IndSym_f,NumberLevelsC,LevelC,
     2        V2Lin_f,KC_f,NumberCoeff_f,tau_l,
     3        J2_i,Eta_i,K_i,V2_i,V1_i,V3_i,NVib_i,Gamma_i,tau_i,
     4        IndSym_i,NumberCoeff_i,NumberLevelsAB,i0,N_i,
     4        V2Lin_i,KC_i,LevelAB,IERROR,theta,NMAX,N_f,N_AB,NMIN,N2XS

      REAL*8  LineStrength,Nu_i_f,Intensity,DeltaTime,G_ns,Minus,
     1        J_i,J_f,Nu0


      INTEGER LeftTime,J2MIN,J2MAX,id

      INTEGER IntJ_i,IntJ_f

      REAL*8 ThreeJSym,SixJSym,EnerSource,A_coef_s_1,A_einst

      INTEGER Date1,Hour1,Min1,Sec1,Date2,Hour2,Min2,Sec2,alloc
      !
      integer StatesUnit,TransUnit
      !
      integer,allocatable :: id_C(:),id_AB(:)
      !
c      INTEGER IndexTerm_AB,IndexTerm_C

C######################################################################C
C#################      START SUBROUTINE     ##########################C
C#################    S(f<--i) = S(C<--AB)   ##########################C

      WRITE(6,3001)
3001  FORMAT(1H0,' ABTOC.INT.INF   START INTENSITY CALCULATIONS '/)


      TempCO=-Planck*VellGT/(Boltz*AbsTemper)
      A_coef_s_1=64.0d-36 * pi**4  / (3.0d0 * planck)

C*********************************************************************
C
C     OPEN FILE FOR STORING UNCONTRACTED COEF. OF A-B
C
C*********************************************************************
       OPEN(UNIT=49,STATUS='SCRATCH',FORM='UNFORMATTED')
C**********************************************************************
C
C     OPEN FILE FOR STORING UNCONTRACTED COEF. OF C STATE
C
C**********************************************************************
      DO i0 = 1,(KMAXC+1)*4
      OPEN(UNIT=100+i0,STATUS='SCRATCH',FORM='UNFORMATTED')
      ENDDO

C**************************************************************************C
C                                                                          C
C     CALCULATES MATRIX ELEMENTS OF THE irreducible repres. D(1)           C
C     OF THE 3D ROT. GROUR                                                 C
C     ON THE ROTATION FUNCTIONS <J K M tau |D_{ss'}^{1} | J' K' M' tau'>   C
C                                                                          C
C**************************************************************************C
C                                                                          C
      WRITE(6,3003)
3003  FORMAT(1H0,'                 CALCUL.  ROT. MATR. ELEMENTS '/)
      CALL Calc_RotateMatrElement(ThreeJsymbol,SixJsymbol)
C                                                                          C
C**************************************************************************C


C**************************************************************************C
C                                                                          C
C     PREREAD EIGENVECTORS UNCONTRACTED COEFFICIENTS
C     FOR LOWER A AND B  STATES AND WRITE TO THE TEMP. FILE 49
C                                                                          C
C**************************************************************************C
C                                                                          C
      WRITE(6,3004)
3004  FORMAT(1H0,'                 PREREAD EIGENVECTORS UNCONTRACTED ',
     1                             'COEFFICIENTS'/)

      CALL PreReadCoeficientsAB(NotateLevelAB,
     1                          EnergyAB,NumberLevelsAB,Q_AB)
C                                                                          C
C**************************************************************************C

C**************************************************************************C
C                                                                          C
C     PREREAD EIGENVECTORS UNCONTRACTED COEFFICIENTS
C     FOR UPPER C STATE AND WRITE TO THE TEMP. FILES 50+J_C
C                                                                          C
C**************************************************************************C
C                                                                          C
      CALL PreReadCoeficientsC (NotateLevelC ,
     1                          EnergyC ,NumberLevelsC ,Q_C )
C                                                                          C
C**************************************************************************C


C**************************************************************************C
C                                                                          C
C     WHAT TIME IS IT NOW?
      CALL CLOCKV
      WRITE (6,5301) NumberLevelsAB
c     WRITE (OutUnit,5300)
C                                                                          C
C**************************************************************************C



        IF ( Te_C+ZeroEnergyC .GT. Te_A+ZeroEnergyA ) THEN   ! ### T0_C > T0_AB ###
          WRITE (TMPUnit,1019)
c    5    Nu_i_f,
c    3    LineStrength,Intensity,A_einst,
c    1    N_f,K_f,KC_f,
c    1    J2_i,N_i,K_i,KC_i,
c    3    GammaLevel(IndSym_f),GammaLevel(IndSym_i),
c    3    V1_f,V2_f,V3_f,V2Lin_f,NVib_f,GammaLevel(Gamma_f*2-1),
c    4    V1_i,V2_i,V3_i,V2Lin_i,NVib_i,
c    3    GammaLevel(Gamma_i*2-1),Eta_i,
c    3    EnergyC(LevelC)+Te_C+ZeroEnergyC,
c    3    EnergyAB((LevelAB))+Te_A+ZeroEnergyA

        WRITE (TMPUnit,1019)


        ELSE !**********                                     ! ### T0_C < T0_AB ###
     
        WRITE (TMPUnit,1019)
        ENDIF


1019  FORMAT("   Nu_i_f","       Sif","          Iif","            Aif",
     c       "           ","J2 Nf Kf KciNi Ki KCGf Gi",
     c       "  V1f V2fV3fVlfNvGv EtV1iV2iV3iVliNvGvi",
     c       "     Ef           Ei")



C**************************************************************************C
C    States file  
C**************************************************************************C

      allocate(id_C(NumberLevelsC),stat = alloc)
      if (alloc/=0) then
        write (out,"(' intens: ',i9,' allocating array -id_C')") alloc
        stop 'id_C - out of memory'
      end if

      allocate(id_AB(NumberLevelsAB),stat = alloc)
      if (alloc/=0) then
        write (out,"(' intens: ',i9,' allocating array -id_AB')") alloc
        stop 'id_AB - out of memory'
      end if


      id = 0
      StatesUnit = 21
      OPEN(UNIT=StatesUnit,
     c     FILE=DirInput(1:NumDirLetters)//
     c     FileOutput(1:NumFileLetters)//".states")

      TransUnit = 22
      OPEN(UNIT=TransUnit,
     c     FILE=DirInput(1:NumDirLetters)//
     c     FileOutput(1:NumFileLetters)//".trans")

     
C**************************************************************************C
C     START LOOPING OVER THE ALL DIFFERENT C-State ENERGY LEVELS
C**************************************************************************C
      DO LevelC = 1,NumberLevelsC            ! ---  Current C-energy level

C**************************************************************************C
C     QUANTUM NUMBERS OF CUREENT A-B ENERGY LEVEL (i)
C**************************************************************************C
C                                                                          C
      K_f    = NotateLevelC(LevelC,2)
      V2_f   = NotateLevelC(LevelC,3)
      V1_f   = NotateLevelC(LevelC,4)
      V3_f   = NotateLevelC(LevelC,5)
      NVib_f = NotateLevelC(LevelC,6)
      Gamma_f= NotateLevelC(LevelC,7)
      J2_f   = NotateLevelC(LevelC,8)
      IndSym_f= NotateLevelC(LevelC,9)        ! --- Symmetry
      NumberCoeff_f= NotateLevelC(LevelC,10)
      N_f    = NotateLevelC(LevelC,11)

      tau_f = MOD(IndSym_f+1,2)

      KC_f = N_f - K_f
      IF ( MOD(tau_f+KC_f,2).NE.0 ) KC_f=KC_f+1
      IF ( KC_f.GT.N_f ) KC_f=KC_f-2
      Minus = -1.0D+00
      IF ( Eta_i.EQ.2 ) Minus = 1.0D+00
      V2Lin_f = 2*V2_f+IABS(K_f+INT(Minus*LambdaC))

          id = id + 1
          id_C(LevelC) = id
          J_f = 0.5D+00*DFLOAT(J2_f)
          
          WRITE (StatesUnit,1021)
     1    id,EnergyC(LevelC)+Te_C,Gns(IndSym_f)*(J2_f+1)*int(XMULTI),
     1    J_f,N_f,K_f,KC_f,
     3    GammaLevel(IndSym_f),
     3    V1_f,V2_f,V3_f,V2Lin_f,NVib_f,GammaLevel(Gamma_f*2-1),3


      ENDDO


1021  FORMAT( i12,1x,f12.6,1x,i6,1x,f7.1,1x,i3,2(1x,i3),1x,
     1        a2,1x,5(1x,i3),1x,a2,1x,i2," ::")



C**************************************************************************C
C                                                                          C
C     START LOOPING OVER THE ALL DIFFERENT A-B ENERGY LEVELS
C                                                                          C
C**************************************************************************C
C                                                                          C
      DO LevelAB = 1,NumberLevelsAB            ! ---  Current AB-energy level


C**************************************************************************C
C                                                                          C
C     QUANTUM NUMBERS OF CUREENT A-B ENERGY LEVEL (i)
C                                                                          C
C**************************************************************************C
C                                                                          C
      Eta_i  = NotateLevelAB(LevelAB,1)
      K_i    = NotateLevelAB(LevelAB,2)
      V2_i   = NotateLevelAB(LevelAB,3)
      V1_i   = NotateLevelAB(LevelAB,4)
      V3_i   = NotateLevelAB(LevelAB,5)
      NVib_i = NotateLevelAB(LevelAB,6)
      Gamma_i= NotateLevelAB(LevelAB,7)
      J2_i   = NotateLevelAB(LevelAB,8)
      IndSym_i= NotateLevelAB(LevelAB,9)        ! --- Symmetry
      NumberCoeff_i= NotateLevelAB(LevelAB,10)
      N_i    = NotateLevelAB(LevelAB,11)

      tau_i = MOD(IndSym_i+1,2)
      KC_i = N_i-K_i
      IF ( MOD(tau_i+KC_i+Eta_i+1,2).NE.0 ) KC_i=KC_i+1
      IF ( KC_i.GT.N_i ) KC_i=KC_i-2
      Minus = -1.0D+00
      IF ( Eta_i.EQ.2 ) Minus = 1.0D+00
      V2Lin_i = 2*V2_i+IABS(K_i+INT(Minus*LambdaAB))

          id_AB(LevelAB) = id
      
          id = id + 1
          J_i = 0.5D+00*DFLOAT(J2_i)
          
          WRITE (StatesUnit,1022)
     1    id,EnergyAB((LevelAB))+Te_A+ZeroEnergyA-ZeroEnergyC,
     1    Gns(IndSym_i)*(J2_i+1),
     1    J_i,N_i,K_i,KC_i,
     3    GammaLevel(IndSym_i),
     4    V1_i,V2_i,V3_i,V2Lin_i,NVib_i,GammaLevel(Gamma_i*2-1),Eta_i
      
      
      ENDDO 

1022  FORMAT( i12,1x,f12.6,1x,i6,1x,f7.1,3(1x,i3),1x,
     1        a2,1x,5(1x,i3),1x,a2,1x,i2" ::")



       close(StatesUnit)
C**************************************************************************C
C                                                                          C
C     START LOOPING OVER THE ALL DIFFERENT A-B ENERGY LEVELS
C                                                                          C
C**************************************************************************C
C                                                                          C
      DO LevelAB = 1,NumberLevelsAB            ! ---  Current AB-energy level


C**************************************************************************C
C                                                                          C
C     QUANTUM NUMBERS OF CUREENT A-B ENERGY LEVEL (i)
C                                                                          C
C**************************************************************************C
C                                                                          C
      Eta_i  = NotateLevelAB(LevelAB,1)
      K_i    = NotateLevelAB(LevelAB,2)
      V2_i   = NotateLevelAB(LevelAB,3)
      V1_i   = NotateLevelAB(LevelAB,4)
      V3_i   = NotateLevelAB(LevelAB,5)
      NVib_i = NotateLevelAB(LevelAB,6)
      Gamma_i= NotateLevelAB(LevelAB,7)
      J2_i   = NotateLevelAB(LevelAB,8)
      IndSym_i= NotateLevelAB(LevelAB,9)        ! --- Symmetry
      NumberCoeff_i= NotateLevelAB(LevelAB,10)
      N_i    = NotateLevelAB(LevelAB,11)

      tau_i = MOD(IndSym_i+1,2)
      KC_i = N_i-K_i
      IF ( MOD(tau_i+KC_i+Eta_i+1,2).NE.0 ) KC_i=KC_i+1
      IF ( KC_i.GT.N_i ) KC_i=KC_i-2
      Minus = -1.0D+00
      IF ( Eta_i.EQ.2 ) Minus = 1.0D+00
      V2Lin_i = 2*V2_i+IABS(K_i+INT(Minus*LambdaAB))
C                                                                          C
C**************************************************************************C


C**************************************************************************C
C                                                                          C
C     READ eigenvector coefficients from temp-file (49)
C                                                                          C
C**************************************************************************C
C                                                                          C
C      CoeffAB(0) = 0.0
      DO TERM = 0,MaxNumCoeffAB
            CoeffAB(term) = 0.0
      ENDDO

      READ(49)( CoeffAB(term), Term = 1,NumberCoeff_i )
C                                                                          C
C**************************************************************************C

C**************************************************************************C
C                                                                          C
C     IndexAB determines correspondens between different indexes
C     eta,K,V2,V1,V3,Nvib,Gamma and so on.
C     and term of CoeffAB matrix
C                                                                          C
C**************************************************************************C
C                                                                          C
      Call DetermCorrespForCoeffAB(J2_i,IndSym_i,NumberCoeff_i,
     1     IndexAB)
C                                                                          C
C**************************************************************************C

C**************************************************************************C
C                                                                          C
C     REWIND ALL FILES WITH UPPER STATE C EIGENFUNCTION
C                                                                          C
C**************************************************************************C
C                                                                          C
      DO i0 = 1,4*(KMAXC+1)
         REWIND(100+i0)
      ENDDO
C                                                                          C
C**************************************************************************C





C**************************************************************************C
C                                                                          C
C     START LOOPING OVER THE ALL DIFFERENT C-State ENERGY LEVELS
C                                                                          C
C**************************************************************************C
C                                                                          C
      DO LevelC = 1,NumberLevelsC            ! ---  Current C-energy level

C**************************************************************************C
C                                                                          C
C     QUANTUM NUMBERS OF CUREENT A-B ENERGY LEVEL (i)
C                                                                          C
C**************************************************************************C
C                                                                          C
      K_f    = NotateLevelC(LevelC,2)
      V2_f   = NotateLevelC(LevelC,3)
      V1_f   = NotateLevelC(LevelC,4)
      V3_f   = NotateLevelC(LevelC,5)
      NVib_f = NotateLevelC(LevelC,6)
      Gamma_f= NotateLevelC(LevelC,7)
      J2_f   = NotateLevelC(LevelC,8)
      IndSym_f= NotateLevelC(LevelC,9)        ! --- Symmetry
      NumberCoeff_f= NotateLevelC(LevelC,10)
      N_f    = NotateLevelC(LevelC,11)

      tau_f = MOD(IndSym_f+1,2)

      KC_f = N_f - K_f
      IF ( MOD(tau_f+KC_f,2).NE.0 ) KC_f=KC_f+1
      IF ( KC_f.GT.N_f ) KC_f=KC_f-2
      Minus = -1.0D+00
      IF ( Eta_i.EQ.2 ) Minus = 1.0D+00
      V2Lin_f = 2*V2_f+IABS(K_f+INT(Minus*LambdaC))

C                                                                          C
C**************************************************************************C


C**************************************************************************C
C                                                                          C
C     QUANTUM NUMBERS OF CUREENT C ENERGY LEVEL (f)
C                                                                          C
C**************************************************************************C
C                                                                          C
C         TAKE FILE WITH UPPER STATE C EIGENFUNCTION
      UNIT_f=100+N_f*4+IndSym_f
C          CoeffC(0) = 0.0
      DO TERM = 0,MaxNumCoeffC
            CoeffC(term) = 0.0
      ENDDO

      READ(Unit_f)( CoeffC(term), Term = 1,NumberCoeff_f )
C                                                                          C
C**************************************************************************C


C**************************************************************************C
C                                                                          C
C     FIRST SELECTION RULES (SYMMETRY)
C                                                                          C
C**************************************************************************C
C                                                                          C
      IF ( ( ( IndSym_i*IndSym_f.EQ.2 .OR.
     1         IndSym_i*IndSym_f.EQ.12 ).AND. SYMM)
     1 .OR. (IndSym_i*IndSym_f.EQ.2 .AND. .NOT.SYMM) )THEN  ! -- A1 <- A2 or B1 <-B2 OR A' <-> A"


C**************************************************************************C
C                                                                          C
C     Transitions under consideration
C                                                                          C
C**************************************************************************C
C                                                                          C
      Nu_i_f = DABS(Te_C-Te_A+ZeroEnergyC-ZeroEnergyA+
     1              EnergyC(LevelC)-EnergyAB(LevelAB) )
      Nu0    = DABS(Te_C+ZeroEnergyC-Te_A-ZeroEnergyA)
C
      IF ( Nu_i_f.LE.NuMax .AND.
     1     Nu_i_f.GE.NuMin ) THEN  ! -- NuMin < Nu < NuMax
C                                                                          C
C**************************************************************************C
C                                                                          C

C**************************************************************************C
C                                                                          C
C     IndexAB determines correspondens between different indexes
C     eta,K,V2,V1,V3,Nvib,Gamma and so on.
C     and term of CoeffC matrix
C                                                                          C
C**************************************************************************C
C                                                                          C
      Call DetermCorrespForCoeffC(N_f,IndSym_f,NumberCoeff_f,IndexC)
C                                                                          C
C**************************************************************************C


C**************************************************************************C
C                                                                          C
C     START CALCULATING OF THE LineStrength VALUE
C                                                                          C
C**************************************************************************C
C                                                                          C
      N2XS=NINT(XMULTI)-1
      NMIN=INT(IABS(J2_i-N2XS)/2)
      NMAX=INT(    (J2_i+N2XS)/2)

      J2MIN=IABS(2*N_f-N2XS)
      J2MAX=     2*N_f+N2XS

      LineStrength = 0.0D+00
      A_einst = 0

      DO J2_f = J2MIN,J2MAX
      IF ( J2_f.GE.IABS(J2_i-2) .AND. J2_f.LE.J2_i+2  ) THEN   ! ---  J_C >= |J_AB-1|

      STemp = 0                           ! -- variable for calculation of the linestrength
      DO N_AB = NMIN,NMAX
      IF ( N_f.GE.IABS(N_AB-1) .AND. N_f.LE.N_AB+1  ) THEN   ! ---  N_f >= |N_AB-1|


C****************************************************************************C
C                                                                            C
       J_i = 0.5D+00*DFLOAT(J2_i)
       J_f = 0.5D+00*DFLOAT(J2_f)
       IntJ_i = INT(J_i)
       IntJ_f = INT(J_f)
       SixJSym = SixJsymbol(IntJ_i+1,N_AB+1,N_f-N_AB+2,IntJ_f-IntJ_i+2)

C                                                                            C
C****************************************************************************C
C                                                                            C
      IF ( SixJSym.NE.0.0 ) THEN

C****************************************************************************C
C                                                                            C
C     START LOOPING OVER THE ALL DIFFERENT NON-ZERO MATRIX ELEMENTS
C                                                                            C
C****************************************************************************C
C                                                                            C
      DO theta =0,1                                          ! --  theta
      DO Term = MinTermK(theta+1),MaxTermK(N_AB+1,theta+1)    ! --  Term
c      DO Term = 1,MaxTermK(KMAXAB+1,2)    ! --  Term


C****************************************************************************C
C                                                                            C
C     DETERMINE OF THE INDEXES FOR THE AB STATE UNCONTRACTING COEFFICIENTS   C
C                                                                            C
C****************************************************************************C
C                                                                            C

      etaAB=  INT(IndexME(1,Term)/10000000)
      Ntemp = IndexME(1,Term)-etaAB*10000000
      K_AB=   INT(Ntemp/100000)
      Ntemp = Ntemp-K_AB*100000
      V2AB  = INT(NTemp/1000)
      Ntemp = Ntemp-V2AB*1000
      NVibAB= INT(NTemp/10)
      GammaVibAB = Ntemp-NVibAB*10

      K_C=INT(IndexME(2,Term)/100000)
      Ntemp = IndexME(2,Term)-K_C*100000
      V2C   = INT(NTemp/1000)
      Ntemp = Ntemp-V2C*1000
      NVibC = INT(NTemp/10)
      GammaVibC = Ntemp-NVibC*10

C                                                                            C
C****************************************************************************C


      IF ( IABS(K_C-K_AB).LE.1
     1    .AND. K_C.LE.N_f .AND. K_AB.LE.N_AB ) THEN  ! -- |K_C-K_AB| <  or = 1


C****************************************************************************C
C                                                                            C
       Minus = -1.0D+00
       IF (Mod(K_AB+N_f+N_AB,2).EQ.0) Minus = 1.0D+00

         ThreeJSym =
     1   ThreeJsymbol( N_AB+1 , N_f-N_AB+2 , K_AB+1 , K_C-K_AB+2 )

C                                                                            C
C****************************************************************************C

         RotMatrElem = Minus*SixJSym*ThreeJSym

C****************************************************************************C
C                                                                            C
C     IF RotMatrElem = 0.0 THEN IT IS USELESS CALCULATE OTHER TERMS
C                                                                            C
C****************************************************************************C
C                                                                            C
      IF ( DABS(RotMatrElem).GT.AlmostZero ) THEN

         MatrSigma =  RotMatrElem*
     1                VibrationMatrixElement(Term)

C****************************************************************************C
C                                                                            C
C     SUMMATION OVER ALL INDEXES FOR CALCULATION OF THE LINESTRENGTHES
C                                                                            C
C****************************************************************************C
C                                                                            C
         TermAB = IndexAB(etaAB,N_AB+1,K_AB+1,V2AB+1,NVibAB,GammaVibAB)
         TermC  = IndexC (             K_C +1,V2C +1,NVibC ,GammaVibC )

      IF ( TermAB.NE.0 .AND. TermC.NE.0 )  THEN
        STemp = STemp + CoeffAB(TermAB)*CoeffC(TermC)*MatrSigma
      ENDIF

C                                                                            C
C****************************************************************************C

      ENDIF           !  --- RotMatrElem <> 0.0
      ENDIF           ! -- |K_C-K_AB| <  or = 1
      ENDDO           ! --- Term
      ENDDO           ! --- theta
      ENDIF           ! --- SixJsimbol
      ENDIF           ! --- |N_f-N_AB| <  or = 1
      ENDDO           ! --- N_ab

C****************************************************************************C
C                                                                            C
C        Linestrength
C                                                                            C
C****************************************************************************C
C                                                                            C
      LineStrength = LineStrength+ 0.25E+00*
     1               (DFLOAT(J2_i)+1.0D+00)*(DFLOAT(J2_f)+1.0D+00)*
     1                STemp*STemp
C                                                                            C
C****************************************************************************C

      ENDIF           ! ---  J_C >= |J_AB-1|
      ENDDO           ! --- J_f
      !
      A_einst = A_coef_s_1*LineStrength*abs(Nu_i_f)**3/
     C          (DFLOAT(J2_f)+1.0D+00)
      !

C****************************************************************************C
C                                                                            C
C        INTENSITY
C                                                                            C
C****************************************************************************C
C                                                                            C

C***********************
        IF ( Te_C+ZeroEnergyC .GT. Te_A+ZeroEnergyA ) THEN   ! ### T0_C > T0_AB ###
             Q = Q_AB+Q_C*DEXP( Nu0*TEMPCO )
             IF ( IntACT ) THEN    ! ### Absorption ###
                 EnerSource = EnergyAB(LevelAB)
                 G_ns=DFloat(Gns(IndSym_i))
             ELSE                  ! ### Emission ###
                 EnerSource = EnergyC (LevelC )
                 G_ns=DFloat(Gns(IndSym_f))
             ENDIF
        ELSE !**********                                     ! ### T0_C < T0_AB ###
             Q = Q_C+Q_AB*DEXP( Nu0*TEMPCO )
             IF ( IntACT ) THEN    ! ### Absorption ###
                 EnerSource = EnergyC (LevelC )
                 G_ns=DFloat(Gns(IndSym_f))
             ELSE  ! *********    ! ### Emission ###
                 EnerSource = EnergyAB(LevelAB)
                 G_ns=DFloat(Gns(IndSym_i))
             ENDIF
        ENDIF
C***********************
C       IF ( .NOT.SYMM ) G_ns= 1.0D+00
C***********************
        IF ( IntACT ) THEN    ! ### Absorption ###
             AbsInt=8.0E-36*PI*PI*PI*
     1              Nu_i_f*(1.0D+00-DEXP(Nu_i_f*TEMPCO))/
     1              (3.0D+00*Planck*VellGT)
         ELSE                  ! ### Emission ###
             AbsInt=64.0E-36*PI*PI*PI*PI*VellGT*
     1              DEXP( Nu0*TEMPCO )*
     1              Nu_i_f*Nu_i_f*Nu_i_f*Nu_i_f/3.0D+00
         ENDIF
C***********************

        Intensity = AbsInt*G_ns*
     1              DEXP(EnerSource*TEMPCO)*
     2              LineStrength/Q


C***********************


C                                                                            C
C****************************************************************************C



C****************************************************************************C
C                                                                            C
C        OUTPUT INTENSITY RESULT
C                                                                            C
C****************************************************************************C
C                                                                            C

        IF ( LineStrength.GE.StrengthLimit .OR.
     1          Intensity.GE.IntensLimit) THEN


        IF ( Te_C+ZeroEnergyC .GT. Te_A+ZeroEnergyA ) THEN   ! ### T0_C > T0_AB ###
          WRITE (TMPUnit,1014)
     5    Nu_i_f,
     3    LineStrength,Intensity,A_einst,
     1    N_f,K_f,KC_f,
     1    J2_i,N_i,K_i,KC_i,
     3    GammaLevel(IndSym_f),GammaLevel(IndSym_i),
     3    V1_f,V2_f,V3_f,V2Lin_f,NVib_f,GammaLevel(Gamma_f*2-1),
     4    V1_i,V2_i,V3_i,V2Lin_i,NVib_i,
     3    GammaLevel(Gamma_i*2-1),Eta_i,
     3    EnergyC(LevelC)+Te_C+ZeroEnergyC-ZeroEnergyA,
     3    EnergyAB((LevelAB))+Te_A

        WRITE (TransUnit,"( i12,1x,i12,1x,1x,es16.8,1x,f16.6,' ||')")
     1         id_C(LevelC),id_AB(LevelAB),A_einst



        ELSE !**********                                     ! ### T0_C < T0_AB ###
          WRITE (TMPUnit,1014)
     5    Nu_i_f,
     3    LineStrength,Intensity,A_einst,
     1    J2_i,N_i,K_i,KC_i,
     1    N_f,K_f,KC_f,
     3    GammaLevel(IndSym_i),GammaLevel(IndSym_f),
     4    V1_i,V2_i,V3_i,V2Lin_i,NVib_i,GammaLevel(Gamma_i*2-1),
     5    Eta_i,
     3    V1_f,V2_f,V3_f,V2Lin_f,NVib_f,GammaLevel(Gamma_f*2-1),
     3    EnergyAB((LevelAB))+Te_A+ZeroEnergyA-ZeroEnergyC,
     3    EnergyC(LevelC)+Te_C

        WRITE (TransUnit,"( i12,1x,i12,1x,1x,es16.8,1x,f16.6,' ||')")
     1         id_AB(LevelAB),id_C(LevelC),A_einst


        ENDIF

1014  FORMAT( 1x,F11.4,1x,D12.4,1x,D15.8,1x,D15.8,
     1        I2,1x,I2,1x,I2,1x,I2,1xI2,1x,I2,1x,I2,1x,
     3        A2,1x,A2,2x,
     1        I2,1x,I2,1x,I2,1x,I2,1x,I2,1x,A2,1x,I1,
     1        I2,1x,I2,1x,I2,1x,I2,1x,I2,1x,A2,1x,
     1        2x,F11.4,2x,F11.4)



c         WRITE (OutUnit,5400)
c    1    J2_f,GammaLevel(IndSym_f),K_f,V2_f,
c    3    GammaLevel(Gamma_f*2-1),NVib_f,EnergyC(LevelC),
c    1    J2_i,GammaLevel(IndSym_i),K_i,V2_i,
c    3    GammaLevel(Gamma_i*2-1),NVib_i,EnergyAB((LevelAB)),
c    5    Nu_i_f,
c    3    LineStrength,Intensity


C                                                                          C
C**************************************************************************C
      ENDIF          ! --- LineStrength <> 0
      ENDIF       ! -- NuMin < Nu < NuMax
      ENDIF       ! -- A1 <- A2 or B1 <-B2
      ENDDO       ! ---  Current C-energy level


C**************************************************************************C
C                                                                          C
C     WHAT TIME IS IT NOW?
      CALL CLOCKV
      DeltaTime  =  DFLOAT(TimePassed)/DFLOAT(LevelAB)
      LeftTime = INT(DeltaTime*NumberLevelsAB)-TimePassed

      Date1 = INT(TimePassed/86400)
      NTemp = TimePassed - Date1*86400
      Hour1 = INT(NTemp/3600)
      NTemp = NTemp - Hour1*3600
      Min1  = INT(NTemp/60)
      Sec1 = NTemp - Min1*60

      Date2 = INT(LeftTime/86400)
      NTemp = LeftTime - Date2*86400
      Hour2 = INT(NTemp/3600)
      NTemp = NTemp - Hour2*3600
      Min2  = INT(NTemp/60)
      Sec2 = NTemp - Min2*60

      WRITE(6,5306) LevelAB,Date1,Hour1,Min1,Sec1,
     1              NumberLevelsAB-LevelAB,Date2,Hour2,Min2,Sec2
C                                                                          C
C**************************************************************************C

      ENDDO            ! ---  Current AB-energy level



C**************************************************************************C
C                                                                          C
C     WHAT TIME IS IT NOW?
      CALL CLOCKV
C                                                                          C
C**************************************************************************C

      DO i0 = 1,4*(2*KMAXC+1)
         CLOSE(100+i0)
      ENDDO
      CLOSE(49)
      CLOSE(TMPUnit)



C****************************************************************************C
C                                                                            C
C        FORMATS
C                                                                            C
C****************************************************************************C
C                                                                            C

1000  FORMAT(A100)
1010  FORMAT(20X,'LINESTRENGTHS FOR C <-- AB '/
     1       '   J`K`G`      E`         J"K"G"      E"',
     2       '    |    S(f<-i)   | ( V1 V2 V3)`N`G` ',
     3       '   eta( V1 V2 V3)"N"G"',' |   nu(f<-i)')

1011  FORMAT( 1X,' ',I2,'',I2,' ',A2,F11.4,'<----',
     1               I2,'',I2,' ',A2,F11.4,' | ',
     4     D12.3,' | (',I2,' ',I2,' ',I2,' )',I2,' ',A2,
     2    '<---- ',I1,'(',I2,' ',I2,' ',I2,' )',I2,' ',A2,' |',F11.4)



1015  FORMAT( 1X,' ',I2,' ',A2,' ',I2,'',I2,'  ',
     1               A2,' ',I2,' ',F11.4,'  ',
     1               I2,'',I2,' ',A2,F11.4/
     1           ' ',I2,' ',A2,' ',I2,'',I2,'  ',
     1               A2,' ',I2,' ',F11.4,'  ',
     4               F11.4,'  ',D12.3,'  ',D12.3)


5300  FORMAT('1',5X,12('*'),' CALCULATED TRANSITIONS ',12('*')//
     1            '0',18X,'J',1X,'SYM.',1X,'KV',1X,'V2',1X,
     1            'V.SYM.',3X,
     1      'N',10X,'ENERGY',8X,'NUCALC',6X,'LINE STRENGTH',
     1      3X,'INTENSITY'//'0',54X,'-1',12X,'-1',16X,'2',
     1      10X,'-1'/
     1           ' ',52X,'CM',12X,'CM',13X,'DEBYE',8X,
     1      'MOL  CM'//' ',128('-')/)


5400  FORMAT('0 UPPER STATE: ',3X,I2,1X,A4,1X,I2,1X,I2,2X,A4,
     1      2X,I2,2X,F15.5/
     1       '  LOWER STATE: ',3X,I2,1X,A4,1X,I2,1X,I2,2X,A4,
     1          2X,I2,2X,2F15.5,2E15.5)

5301  FORMAT(1H0,' ABTOC.INT.INF   '
     1           ' UNDER CONSIDIRATION ',I8,'  AB Energy levels ')

5305  FORMAT('1     ',A8)

5306  FORMAT(1H0,'                 ',
     1       I5,' (',I2,':',I2,':',I2,':',I2,')  ---> ',
     1       I5,' (',I2,':',I2,':',I2,':',I2,')' )


C                                                                            C
C****************************************************************************C
C                                                                            C



      RETURN


1009      WRITE(*,*) 'ERROR: FILE ',IERROR,' WAS NOT OPENED'

          STOP



      END
C####################       Calc_LineStrength       ########################
C###########################################################################



C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : READ EIGENVECTORS UNCONTRACTED COEFFICIENTS
C             : FOR LOWER A and B STATES
C             :
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################

      SUBROUTINE DetermCorrespForCoeffAB(J2_AB,IndSym,
     1           NumberCoeffAB,IndexAB)


C     IMPLICIT REAL*8 (A-H,O-Z)
C     IMPLICIT LOGICAL (A-Z)

      CHARACTER*200 LINE
      CHARACTER*2 TOTSYM

      include 'mu.h'
      include 'dipolsys.h'
      include 'molcul.h'

      INTEGER K,ISurf,GammaVib,NVib0,V2,Term,NMIN,NMax,
     1        NumberCoeffAB,tau,Flag,IndSym,J2_AB,N2XS,N


      INTEGER IndexAB(2,KMAXAB+1,KMAXAB+1,NV2MAXAB+1,MaxNvibMax,2)



C ---------------------------------------------------------- C
C --------------------  A and B states --------------------- C
C                                                                                       C

      IF (SYMM) THEN
        IF (IndSym .EQ. 1) THEN
           tau=0
           FLAG=0
        ENDIF
        IF (IndSym .EQ. 2) THEN
           tau=1
           FLAG=1
        ENDIF
        IF (IndSym .EQ. 4) THEN
           tau=1
           FLAG=0
        ENDIF
        IF (IndSym .EQ. 3) THEN
           tau=0
           FLAG=1
        ENDIF
      ELSE
        IF (IndSym .EQ. 1) THEN
           tau=0
           FLAG=0
        ENDIF
        IF (IndSym .EQ. 2) THEN
           tau=1
           FLAG=0
        ENDIF
      ENDIF

      N2XS=XMULTI-1

      NMIN=IABS(INT( (J2_AB-N2XS)/2 ) )
      NMAX=     INT( (J2_AB+N2XS)/2 )


C     --- MAKE IndexAB TO BE A ZERO

      DO N=0,KMAXAB
      DO K=0,KMAXAB
      DO ISurf=1,2,1
      DO V2=0,NV2MAX(ISurf)
      DO GammaVib=1,2
      DO NVib0 = 1,MaxNvibMax,1  ! -- NVib0
           IndexAB(ISurf,N+1,K+1,V2+1,NVib0,GammaVib) = 0
      ENDDO                ! --- NVib0
      ENDDO                ! --- GammaVib
      ENDDO                ! --- V2
      ENDDO  ! ISurf
      ENDDO  ! k
      ENDDO  ! N

      Term  = 0       ! --  Term = 0..NumberCoeffAB
C
C     LOOP OVER DIFFERENT K QUANTUM NUMBERS AND LOWER STATES
C

      DO N=NMIN,NMAX
      DO K=0,N
      DO ISurf=1,2,1
C
C        NOTE: N+tau always even for K=0
C

      IF ( K .EQ.0 .AND. MOD(tau+N+ISurf+1,2).EQ.0 .OR.
     1     K .NE.0) THEN
C
C        Determine symmetry of the stretching function
C
         IF (MOD(FLAG+K,2).EQ.0) THEN
            GammaVib=1
         ELSE
            GammaVib=2
         ENDIF
         IF ( .NOT.SYMM ) GammaVib=1

         DO V2=0,NV2MAX(ISurf)
           DO NVib0 = 1,Nvib(ISurf,GammaVib),1  ! -- NVib0

           Term = Term +1
           IndexAB(ISurf,N+1,K+1,V2+1,NVib0,GammaVib) = Term

         ENDDO                ! --- NVib0
         ENDDO                ! --- V2
      ENDIF  ! tau+J = even

      ENDDO  ! ISurf
      ENDDO  ! k
      ENDDO  ! N


      IF ( Term .NE. NumberCoeffAB ) THEN
          WRITE(*,*) 'ERROR: INPUT DATA DO NOT',
     1                  ' CORRESPOND TO MAX PARAM-Rs 1'
            STOP

      ENDIF

C                                                            C
C --------------------  A and B states --------------------- C
C ---------------------------------------------------------- C

      RETURN

      END

C#########################      ReadCoeficientsAB     ######################
C###########################################################################


C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : READ EIGENVECTORS UNCONTRACTED COEFFICIENTS
C             : FOR UPPER C STATE
C             :
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################

      SUBROUTINE DetermCorrespForCoeffC(N_C,IndSym,NumberCoeffC,IndexC)

      include 'mu.h'
      include 'dipolsys.h'
      include 'molcul.h'

      CHARACTER*200 LINE
      INTEGER K,GammaVib,NVib0,V2,Term,
     1        NumberCoeffC,tau,Flag,IndSym,
     1        N_C

      INTEGER IndexC(KMAXC+1,NV2MAXC+1,MaxNvibMax,2)


C ---------------------------------------------------------- C
C --------------------         C       --------------------- C


      IF (SYMM) THEN
        IF (IndSym .EQ. 1) THEN
           tau=0
           FLAG=0
        ENDIF
        IF (IndSym .EQ. 2) THEN
           tau=1
           FLAG=1
        ENDIF
        IF (IndSym .EQ. 4) THEN
           tau=1
           FLAG=0
        ENDIF
        IF (IndSym .EQ. 3) THEN
           tau=0
           FLAG=1
        ENDIF
      ELSE
        IF (IndSym .EQ. 1) THEN
           tau=0
           FLAG=0
        ENDIF
        IF (IndSym .EQ. 2) THEN
           tau=1
           FLAG=0
        ENDIF
      ENDIF

C     --- MAKE IndexC TO BE A ZERO

      DO K=0,KMAXC
      DO V2=0,NV2MAXC
      DO GammaVib=1,2
      DO NVib0 = 1,MaxNvibMax,1  ! -- NVib0
           IndexC(K+1,V2+1,NVib0,GammaVib) = 0
      ENDDO                ! --- NVib0
      ENDDO                ! --- GammaVib
      ENDDO                ! --- V2
      ENDDO  ! k


      Term  = 0       ! --  Term = 0..NumberCoeffC

C
C     LOOP OVER DIFFERENT EVEN K QUANTUM NUMBERS
C
      DO K=0,N_C,2     !  k - even
C
C        NOTE: N+tau always even for K=0
C
      IF ( K .EQ.0 .AND. MOD(tau+N_C,2).EQ.0 .OR.
     1     K .NE.0) THEN

C
C        Determine symmetry of the stretching function
C
         IF ( MOD( FLAG,2 ).EQ.0 ) THEN
            GammaVib=1
         ELSE
            GammaVib=2
         ENDIF
         IF (.NOT.SYMM) GammaVib=1

         DO V2=0,NV2MAXC
           DO NVib0 = 1,Nvib(3,GammaVib),1  ! -- NVib0
           Term = Term +1
           IndexC(K+1,V2+1,NVib0,GammaVib) = Term

         ENDDO                ! --- NVib0
         ENDDO                ! --- V2
      ENDIF  ! tau+J = even
      ENDDO  ! k - even

C
C     LOOP OVER DIFFERENT ODD K QUANTUM NUMBERS
C
      DO K=1,N_C,2     !  k - even
C
C        Determine symmetry of the stretching function
C
         IF (MOD(FLAG+1,2).EQ.0) THEN
            GammaVib=1
         ELSE
            GammaVib=2
         ENDIF
         IF (.NOT.SYMM) GammaVib=1

         DO V2=0,NV2MAXC
           DO NVib0 = 1,Nvib(3,GammaVib),1  ! -- NVib0
           Term = Term +1
           IndexC(K+1,V2+1,NVib0,GammaVib) = Term

         ENDDO                ! --- NVib0
         ENDDO                ! --- V2
      ENDDO  ! k - odd

      IF ( Term .NE. NumberCoeffC ) THEN
          WRITE(*,*) 'ERROR: INPUT DATA DO NOT',
     1               ' CORRESPOND TO MAX PARAM-Rs 2',
     1               Term,NumberCoeffC
            STOP

      ENDIF

C                                                            C
C --------------------     C states    --------------------- C
C ---------------------------------------------------------- C



      RETURN

      END

C########################      ReadCoeficientsC       ######################
C###########################################################################




C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : PREREAD EIGENVECTORS UNCONTRACTED COEFFICIENTS
C             : FOR LOWER A and B STATES
C             :
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################
      SUBROUTINE PreReadCoeficientsAB(NotateLevelAB,EnergyAB,
     1                                NumberLevelsAB,Q_AB)


      CHARACTER*24 LINE
      CHARACTER*2 TOTSYM

      include 'mu.h'
      include 'dipolsys.h'
      include 'molcul.h'

      INTEGER K,ISurf,GammaVib,NVib0,V2,Term,
     1        NumberCoeffAB,tau,Flag,IndSym,J2_AB,
     2        IndSymAB,NotateLevelAB(MaxNumLevelAB,11)

      INTEGER MaxK,MaxEta,MaxV2,MaxGamma,MaxNVib,MaxN

      INTEGER Level,NumberLevelsAB,N2XS,N,NMIN,NMAX

      REAL*8 CoeffAB(0:MaxNumCoeffAB),CoeffMax

      REAL*8 EnergyAB(MaxNumLevelAB),Energy,Q_AB,G_ns

      CHARACTER*80 FILNAMAB


C     OPEN FILE WITH LOWER STATE A,B EIGENFUNCTION
C
      FILNAMAB=DirInput(1:NumDirLetters)//FileCoefAB
      OPEN(UNIT=94,FILE=FILNAMAB)


      Level = 0

      READ(94,*) ZeroEnergyA

      Q_AB = 0

      N2XS=INT(XMULTI)-1

      DO WHILE (.NOT. EOF(94))



C
C     SKIP FIRST LINE
C
         READ(94,1000) LINE
C
C     READ NEXT LABELING LINES
C

      READ(94,*) J2_AB
      READ(94,502) TOTSYM
      READ(94,*) Energy


      NMIN=INT(IABS(J2_AB-N2XS)/2)
      NMAX=INT(    (J2_AB+N2XS)/2)



C      READ(94,501) JX2VAL
C      READ(94,502) TOTSYM
C      READ(94,503) EnergyAB

      if (SYMM) then 
        SELECT CASE (TOTSYM)
        case ('A1')
           IndSym = 1
           tau=0
           FLAG=0
        case ('A2')
           IndSym = 2
           tau=1
           FLAG=1
        case ('B1')
           IndSym = 4
           tau=1
           FLAG=0
        case ('B2')
           IndSym = 3
           tau=0
           FLAG=1
        case default
          write(6,'("Illegal TOTSYM")')
          stop 'illegal TOTSYM'
        end select
        !
      else
        !
        SELECT CASE (TOTSYM)
        !
        case ('A''','A1')
           IndSym = 1
           tau=0
           FLAG=0
        case ('A"','B1')
           IndSym = 2
           tau=1
           FLAG=0
        case default
          write(6,'("Illegal TOTSYM")')
          stop 'illegal TOTSYM'
        end select
        !
      endif
      G_ns=DFloat(Gns(IndSym))
      if ( .NOT.SYMM ) G_ns= 1.0D+00

      IF ( (Energy*TEMPCO).GT.-100 ) THEN

      Q_AB = Q_AB + G_ns*DFLOAT(J2_AB+1)*
     1              DEXP(Energy*TEMPCO)

      ENDIF

C
C     read eigenvector coefficients
C

      NumberCoeffAB = 0
      CoeffAB(0) = 0.0d+00
      DO WHILE (CoeffAB(NumberCoeffAB).GE.(-1.0))
            NumberCoeffAB = NumberCoeffAB +1
          READ(94,503) CoeffAB(NumberCoeffAB)
      ENDDO
      NumberCoeffAB = NumberCoeffAB - 1

C---- CHECK THE INPUT DATA
C     --- TRANSITIONS ONLY FROM LEVELS, LOWER THEN EnerMaxAB LIMIT
      IF ( Energy.LE.EnerMaxAB ) THEN

      Term  = 0       ! --  Term = 0..NumberCoeffAB
C
C     LOOP OVER DIFFERENT K QUANTUM NUMBERS AND LOWER STATES
C

C     ------ PARAMETERS FOR FINDING MAX EIGENVECTS COEFF-S
      CoeffMax = DAbs(CoeffAB(0))
      MaxN = 0
      MaxK = 0
      MaxEta = 1
      MaxV2 = 0
      MaxGamma = 1
      MaxNVib = 1

      DO N=NMIN,NMAX,1
      DO K=0,N
      DO ISurf=1,2,1
C
C        NOTE: N+tau always even for K=0
C

      IF ( ( K .EQ.0 .AND. MOD(tau+N+ISurf+1,2).EQ.0 ) .OR.
     1     ( K .NE.0 ) ) THEN
C
C        Determine symmetry of the stretching function
C
         IF (MOD(FLAG+K,2).EQ.0) THEN
            GammaVib=1
         ELSE
            GammaVib=2
         ENDIF
         IF (.NOT.SYMM) GammaVib=1

         DO V2=0,NV2MAX(ISurf)
           DO NVib0 = 1,Nvib(ISurf,GammaVib),1  ! -- NVib0

           Term = Term +1


C     --------  SEARCHING OF THE MAX CoeffAB
           IF ( CoeffMax.LT.DABS(CoeffAB(Term)) ) THEN
               MaxN     = N
               MaxK     = K
               MaxEta   = ISurf
               MaxV2    = V2
               MaxGamma = GammaVib
               MaxNVib  = NVib0
               CoeffMax = DABS(CoeffAB(Term))
           ENDIF
         ENDDO                ! --- NVib0
         ENDDO                ! --- V2
      ENDIF  ! tau+J = even

      ENDDO  ! ISurf
      ENDDO  ! k
      ENDDO  ! N

      IF ( Term .NE. NumberCoeffAB ) THEN
          WRITE(*,*) 'ERROR: INPUT DATA DO NOT',
     1                  ' CORRESPOND TO MAX PARAM-Rs 3'
            STOP

      ENDIF


C     --- TRANSITIONS ONLY FROM THIS INTERVAL OF NVib0 AND GamVib OF A-B PAIR
      IF ( MaxEta .GE.MinSurf      .AND.
     1     MaxEta .LE.MaxSurf      .AND.     !  ---  [Min<=NvibAB<=Max]
     1     MaxNVib.GE.MinNVibAB    .AND.
     1     MaxNVib.LE.MaxNVibAB    .AND.     !  ---  [Min<=NvibAB<=Max]
     2     MaxGamma.GE.MinGamVibAB .AND.
     3     MaxGamma.LE.MaxGamVibAB .AND.     !  ---  [Min<=GamVibAB<=Max]
     3     MaxV2.GE.MinV2AB        .AND.
     3     MaxV2.LE.MaxV2AB        .AND.     !  ---  [Min<=V2AB<=Max]
     3     J2_AB.GE.2*MinJAB       .AND.
     3     J2_AB.LE.2*MaxJAB       .AND.     !  ---  [Min<=JAB<=Max]
     3     MaxK.GE.MinKAB          .AND.
     3     MaxK.LE.MaxKAB      )    THEN     !  ---  [Min<=KAB<=Max]


      Level = Level + 1

      IF ( Level .GT. MaxNumLevelAB ) THEN
        WRITE(*,*) 'ERROR: PARAMETER MaxNumLevelAB',
     1  ' TOO SMALL TO BE A BOUND OF THE  ARRAY',
     1  ' NotateLevelAB, LESS THAN ',Level
          STOP
      ENDIF




      NotateLevelAB(Level,1) = MaxEta
      NotateLevelAB(Level,2) = MaxK
      NotateLevelAB(Level,3) = MaxV2
      NotateLevelAB(Level,4) = INT(NotateStretch
     1                            (MaxEta,MaxGamma,MaxNVib)/100)                 ! --  V1
      NotateLevelAB(Level,5) = NotateStretch(MaxEta,MaxGamma,MaxNVib)-           ! --  V3
     1                         NotateLevelAB(Level,4)*100
      NotateLevelAB(Level,6) = MaxNVib
      NotateLevelAB(Level,7) = MaxGamma
      NotateLevelAB(Level,8) = J2_AB
      NotateLevelAB(Level,9) = IndSym
      NotateLevelAB(Level,10)= NumberCoeffAB
      NotateLevelAB(Level,11) = MaxN

      EnergyAB(Level) = Energy

C            WRITE(51,REC=Level) ( CoeffAB(term), Term = 1,NumberCoeffAB  )
            WRITE(49) ( CoeffAB(term), Term = 1,NumberCoeffAB  )



C******************************************************************************
C*         print out the uncontracted coeffitients                            *
      IF ( OUT_COEFF.EQ.1 )  THEN

           Term = 0
          write(6,*) Energy

          write(6,111) Energy,J2_AB,MaxEta,MaxK,MaxV2,MaxNVib,
     1             GammaLevel(MaxGamma*2-1),
     1             GammaLevel(IndSym)

111   FORMAT(D16.6,2X,I3,1X,I3,1X,I3,1X,I3,I3,2X,A2,2X,A2)


      DO N=NMIN,NMAX,1
      DO K=0,N
      DO ISurf=1,2,1
C
C        NOTE: N+tau always even for K=0
C

      IF ( ( K .EQ.0 .AND. MOD(tau+N+ISurf+1,2).EQ.0 ) .OR.
     1     ( K .NE.0 ) ) THEN
C
C        Determine symmetry of the stretching function
C
         IF (MOD(FLAG+K,2).EQ.0) THEN
            GammaVib=1
         ELSE
            GammaVib=2
         ENDIF
         IF (.NOT.SYMM) GammaVib=1

         DO V2=0,NV2MAX(ISurf)
           DO NVib0 = 1,Nvib(ISurf,GammaVib),1  ! -- NVib0

           Term = Term +1

          IF ( DABS(CoeffAB(Term)).GE.0.0001 ) THEN
          write(6,111) CoeffAB(Term),N,ISurf,K,V2,NVib0,
     1             GammaLevel(GammaVib*2-1),
     1             GammaLevel(IndSym)
          ENDIF

         ENDDO                ! --- NVib0
         ENDDO                ! --- V2
      ENDIF  ! tau+J = even

      ENDDO  ! ISurf
      ENDDO  ! k
      ENDDO  ! N

      ENDIF   ! --- OUT_COEFF
C******************************************************************************


      ENDIF ! ---  Indexes are in interval under consideration
      ENDIF ! ---  Energy<EnerMaxAB
      ENDDO ! ---  FILE AB
      REWIND(49)



      NumberLevelsAB = Level

      CLOSE(94)



C     INTEGER I0,I1
C     REAL*8 tttt(10,10)
C      OPEN(UNIT=1,ACCESS='DIRECT',
C     1     RECL=8*10)
C
C      DO i1 = 1,10
C          WRITE(1,REC=I1) (SQRT(1.0*I1*i0),i0 = 1,10)
C      ENDDO
C
C
CC       CLOSE(1)
C
C
CC      REWIND 1
C
CC       OPEN(UNIT=1,ACCESS='DIRECT',
CC     1     RECL=80)
C
C
C      DO i1 = 1,10
C          READ(1,REC=I1) (TTTT(I1,I0),i0 = 1,10)
C      ENDDO
C
C
C
C      CLOSE(1)




C                                                            C
C --------------------  A and B states --------------------- C
C ---------------------------------------------------------- C


1000  FORMAT(A24)
502   FORMAT(A2)
503   FORMAT(D24.14)


      RETURN

      END

C#######################      PreReadCoeficientsAB     #####################
C###########################################################################



C##########################################################################
C             :
C DATE        : 05.06.2000
C AUTHOR      : Yurchenko
C UPDATES     :
C LANGUAGE    : FORTRAN
C PURPOSE     : PREREAD EIGENVECTORS UNCONTRACTED COEFFICIENTS
C             : FOR UPPER C STATE
C             :
C SUBPROGRAMS :
C      CALLED :
C             :
C###########################################################################

      SUBROUTINE PreReadCoeficientsC(NotateLevelC,EnergyC,
     1                               NumberLevelsC,Q_C)

C     IMPLICIT REAL*8 (A-H,O-Z)
C     IMPLICIT LOGICAL (A-Z)

      include 'mu.h'
      include 'dipolsys.h'
      include 'molcul.h'

      CHARACTER*200 LINE
      INTEGER K,GammaVib,NVib0,V2,Term,
     1        NumberCoeffC,tau,Flag,IndSym,
     1        J2_C,N_C

      INTEGER NotateLevelC(MaxNumLevelC,11),Level,NumberLevelsC

      REAL*8 CoeffC(0:MaxNumCoeffC),CoeffMax,Energy

      REAL*8 EnergyC(MaxNumLevelC),Q_C,G_ns

      INTEGER MaxK,MaxV2,MaxGamma,MaxNVib

      CHARACTER*80 FILNAMC


C ---------------------------------------------------------- C
C --------------------         C       --------------------- C
C                                                                                       C

      FILNAMC=DirInput(1:NumDirLetters)//FileCoefC
      OPEN(UNIT=91,FILE=FILNAMC)

      Level = 0

      READ(91,*) ZeroEnergyC

      Q_C = 0.0


      READ(91,1000) LINE
      DO WHILE (line(1:3)/='eof')

C
C     SKIP FIRST LINE
C
C
C     READ NEXT LABELING LINES
C

      READ(91,*) Energy
      READ(91,*) N_C
      READ(91,*) tau
      READ(91,*) Flag
      READ(91,1000) LINE        ! ---  SKIP NEXT LINE WITH V2MAX
      READ(91,*) NumberCoeffC
! ---  SKIP NEXT 4 LINES WITH IROTE,IROTOU,ICEVEU,ICODDU
      DO term=1,4
         READ(91,1000) LINE
      ENDDO
      !
      !
      J2_C = 2*N_C+int(XMULTI)-1


C     ---     determine the total symmetry

      IF (SYMM) THEN
        IF     (tau .EQ.0  .AND. Flag.EQ.0 ) THEN
                IndSym = 1
        ELSEIF (tau .EQ.1  .AND. Flag.EQ.1 ) THEN
                IndSym = 2
        ELSEIF (tau .EQ.1  .AND. Flag.EQ.0 ) THEN
                IndSym = 4
        ELSE
                IndSym = 3
        ENDIF
      ELSE
        IF     (tau .EQ.0  .AND. Flag.EQ.0 ) THEN
                IndSym = 1
        ELSEIF (tau .EQ.1  .AND. Flag.EQ.0 ) THEN
                IndSym = 2
        ENDIF
      ENDIF
      !
      G_ns=DFloat(Gns(IndSym))
      !
      IF ( (Energy*TEMPCO).GT.-100 ) THEN
        Q_C = Q_C + G_ns*( 2.0D+00*DFLOAT(N_C) +1.0D+00)*XMULTI*
     1              DEXP(Energy*TEMPCO)
      ENDIF

C
C     read eigenvector coefficients
C
      DO Term = 1,NumberCoeffC
          READ(91,503) CoeffC(Term)
      ENDDO


C---- CHECK THE INPUT DATA
C     --- TRANSITIONS ONLY FROM LEVELS, LOWER THEN EnerMaxC LIMIT
      IF ( Energy.LE.EnerMaxC ) THEN

      Term  = 0       ! --  Term = 0..NumberCoeffC

      CoeffMax = DAbs(CoeffC(1))
      MaxK = 0
      MaxV2 = 0
      MaxGamma = 1
      MaxNVib = 1


C
C     LOOP OVER DIFFERENT EVEN K QUANTUM NUMBERS
C
      DO K=0,N_C,2     !  k - even
C
C        NOTE: N+tau always even for K=0
C
      IF ( K .EQ.0 .AND. MOD(tau+N_C,2).EQ.0 .OR.
     1     K .NE.0) THEN

C
C        Determine symmetry of the stretching function
C
         IF ( MOD( FLAG,2 ).EQ.0 ) THEN
            GammaVib=1
         ELSE
            GammaVib=2
         ENDIF
         IF (.NOT.SYMM) GammaVib=1

         DO V2=0,NV2MAXC
           DO NVib0 = 1,Nvib(3,GammaVib),1  ! -- NVib0
           Term = Term +1

           IF ( CoeffMax.LT.DABS(CoeffC(Term)) ) THEN
               MaxK     = K
               MaxV2    = V2
               MaxGamma = GammaVib
               MaxNVib  = NVib0
               CoeffMax = DABS(CoeffC(Term))
           ENDIF

         ENDDO                ! --- NVib0
         ENDDO                ! --- V2
      ENDIF  ! tau+J = even
      ENDDO  ! k - even

C
C     LOOP OVER DIFFERENT ODD K QUANTUM NUMBERS
C
      DO K=1,N_C,2     !  k - even
C
C        Determine symmetry of the stretching function
C
         IF (MOD(FLAG+1,2).EQ.0) THEN
            GammaVib=1
         ELSE
            GammaVib=2
         ENDIF
         IF (.NOT.SYMM) GammaVib=1

         DO V2=0,NV2MAXC
           DO NVib0 = 1,Nvib(3,GammaVib),1  ! -- NVib0
           Term = Term +1

           IF ( CoeffMax.LT.DABS(CoeffC(Term)) ) THEN
               MaxK     = K
               MaxV2    = V2
               MaxGamma = GammaVib
               MaxNVib  = NVib0
               CoeffMax = DABS(CoeffC(Term))
           ENDIF

         ENDDO                ! --- NVib0
         ENDDO                ! --- V2
      ENDDO  ! k - odd




      IF ( Term .NE. NumberCoeffC ) THEN
          WRITE(*,*) 'ERROR: INPUT DATA DO NOT',
     1               ' CORRESPOND TO MAX PARAM-Rs 2',
     1               Term,NumberCoeffC,Energy
            STOP

      ENDIF




C     --- TRANSITIONS ONLY FROM THIS INTERVAL OF NVib0 AND GamVib OF C STATE
      IF ( MaxNVib.GE.MinNVibC    .AND.
     1     MaxNVib.LE.MaxNVibC    .AND.     !  ---  [Min<=NvibC<=Max]
     2     MaxGamma.GE.MinGamVibC .AND.
     3     MaxGamma.LE.MaxGamVibC .AND.     !  ---  [Min<=GamVibC<=Max]
     3     MaxV2.GE.MinV2C        .AND.
     3     MaxV2.LE.MaxV2C        .AND.     !  ---  [Min<=V2C<=Max]
     3     J2_C.GE.2*MinJC        .AND.
     3     J2_C.LE.2*MaxJC        .AND.     !  ---  [Min<=JC<=Max]
     3     MaxK.GE.MinKC          .AND.
     3     MaxK.LE.MaxKC      )    THEN     !  ---  [Min<=KC<=Max]


      Level = Level + 1

      IF ( Level .GT. MaxNumLevelC ) THEN
        WRITE(*,*) 'ERROR: PARAMETER MaxNumLevelC',
     1  ' TOO SMALL TO BE A BOUND OF THE  ARRAY',
     1  ' NotateLevelC, LESS THAN ',Level
          STOP
      ENDIF



      NotateLevelC(Level,1) = 3
      NotateLevelC(Level,2) = MaxK
      NotateLevelC(Level,3) = MaxV2
      NotateLevelC(Level,4) = INT(NotateStretch
     1                            (3,MaxGamma,MaxNVib)/100)                 ! --  V1
      NotateLevelC(Level,5) = NotateStretch(3,MaxGamma,MaxNVib)-            ! --  V3
     1                         NotateLevelC(Level,4)*100
      NotateLevelC(Level,6) = MaxNVib
      NotateLevelC(Level,7) = MaxGamma
      NotateLevelC(Level,8) = J2_C
      NotateLevelC(Level,9) = IndSym
      NotateLevelC(Level,10)= NumberCoeffC
      NotateLevelC(Level,11)= N_C

      EnergyC(Level) = Energy

      WRITE(100+N_C*4+IndSym) ( CoeffC(term), Term = 1,NumberCoeffC  )


C******************************************************************************
C*         print out the uncontracted coeffitients                            *
      IF ( OUT_COEFF.EQ.1 )  THEN


C*****************************************************
C     LOOP OVER DIFFERENT EVEN K QUANTUM NUMBERS
           Term = 0
          write(6,*) Energy

          write(6,111) Energy,N_C,MaxK,MaxV2,MaxNVib,
     1             GammaLevel(MaxGamma*2-1),
     1             GammaLevel(IndSym)

      DO K=0,N_C,2     !  k - even
      IF ( K .EQ.0 .AND. MOD(tau+N_C,2).EQ.0 .OR.
     1     K .NE.0) THEN

         IF ( MOD( FLAG,2 ).EQ.0 ) THEN
            GammaVib=1
         ELSE
            GammaVib=2
         ENDIF
         IF (.NOT.SYMM) GammaVib=1

         DO V2=0,NV2MAXC
           DO NVib0 = 1,Nvib(3,GammaVib),1  ! -- NVib0
           Term = Term +1

          IF ( DABS(CoeffC(Term)).GE.0.0001 ) THEN
          write(6,111) CoeffC(Term),N_C,K,V2,NVib0,
     1             GammaLevel(GammaVib*2-1),
     1             GammaLevel(IndSym)
          ENDIF

         ENDDO                ! --- NVib0
         ENDDO                ! --- V2
      ENDIF  ! tau+J = even
      ENDDO  ! k - even

      DO K=1,N_C,2     !  k - even
         IF (MOD(FLAG+1,2).EQ.0) THEN
            GammaVib=1
         ELSE
            GammaVib=2
         ENDIF
         IF (.NOT.SYMM) GammaVib=1

         DO V2=0,NV2MAXC
           DO NVib0 = 1,Nvib(3,GammaVib),1  ! -- NVib0
           Term = Term +1

          IF ( DABS(CoeffC(Term)).GE.0.0001 ) THEN

          write(6,111) CoeffC(Term),N_C,K,V2,NVib0,
     1             GammaLevel(GammaVib*2-1),
     1             GammaLevel(IndSym)

          ENDIF
         ENDDO                ! --- NVib0
         ENDDO                ! --- V2
      ENDDO  ! k - odd

      ENDIF
C****************************************************************************

      ENDIF ! ---  Indexes are in interval under consideration

      ENDIF ! ---  Energy<EnerMaxC


      READ(91,1000,end=118) LINE
      !
      CYCLE
      !
  118 EXIT

      ENDDO ! ---  FILE C
      CLOSE(91)

      NumberLevelsC = Level


C                                                            C
C --------------------     C states    --------------------- C
C ---------------------------------------------------------- C


1000  FORMAT(A100)


503   FORMAT(D24.14)
111   FORMAT(D16.6,2X,I3,'-',I3,'-',I3,I3,2X,A2,2X,A2)


      RETURN

      END

C######################      PeReadCoeficientsC       ######################
C###########################################################################
