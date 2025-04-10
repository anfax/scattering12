!!!*************************************************************
! 文件/File: datain.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: datain.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE DATAIN
C
C     DATAIN - read in data and set up constants
C
C     Modified 3/18/91 to include centrifugal oscillator bend energies
C
C  Called by:
C     VTSTM - main program
C
C  Calls:
C     PREPEF - set up potential parameters
C     PEF    - evaluate potential
C     BEND   - compute bending potential parameters
C     D2VDU2 - compute second derivative along u coordinate
C     ELPART - compute electronic partition functions
C     GAUSSQ - evaluate gaussian quadrature weights and nodes
C     GEOM   - find equilibrium geometry for reactants and products
C     PFBR   - computes classical bend-rotational (coupled) partition
C              functions
C     PFCNB  - compute bending part. fcn.
C     PFCNR  - compute rotational part. fcn.
C     PFCNST - compute stretching part. fcn.
C     SADDLE - find saddle point(s) and do normal mode analysis
C     TITLES - print out 3-line title
C     VMIN   - find minimum energy along u coordinate
C     WKBRD  - read in WKB energies, etc. from unit 14
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSAD
      CHARACTER*9 AA(3)
      CHARACTER*10 AB(2)
      CHARACTER*2 AX(3)
      CHARACTER*5 ADUMX(3)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)
      COMMON /ADIAB1/ VAD(NSDM,2), VMAX(2), VR, VP, SMAX(2), ISMAX(2)
      COMMON /ADIAB2/ SAD1, SAD2, NPAD
      COMMON /BENDOP/ THETA(2), A11, A12, A21, A22, IBOPT
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     * NSSP(NSADDM), LRSTRT
      COMMON /BENDTS/ FBTS(NSADDM), QFBTS(NSADDM), GBTS(NSADDM),
     * XMOMTS(NSADDM)
      COMMON /COBND/  NLEVEL
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /C12/    C12, CONC12, INDFR
      COMMON /CURVE1/ DELCUR, DELCSP, SM, SP
      COMMON /CURVE2/ SM2, SP2
      COMMON /EACT1/   IACT, NT1(10), NT2(10)
      LOGICAL LMAX
      COMMON /EBND1/   LMAX
      PARAMETER (NTEMDM=100)                                            TCA0996
      DOUBLE PRECISION KAPW
      LOGICAL PTEMP                                                     TCA0996
      COMMON /TEMPCM/ TEMP(NTEMDM), BETA(NTEMDM), CPHI(NTEMDM),
     * CNST(NTEMDM), CPHIC(NTEMDM), CNSTC(NTEMDM), EQUIL(NTEMDM,2),
     * KAPW(NTEMDM,2), RATIO(NTEMDM,2), NTMAX, PTEMP(NTEMDM)            TCA0996
      DIMENSION EFACTF(NTEMDM), EFACTR(NTEMDM), ITEMP(NTEMDM)           TCA0996
      COMMON /ELFCT/  ELFACT(NTEMDM,2)
      COMMON /ESPEC/  ESPEC(40), NESPEC
      COMMON /GTST1/   NPGTST
      COMMON /INTER1/ NINT
      LOGICAL LLAG, LLAGRS
      COMMON /LAGCOM/ PT3(NQGKDM), WT3(NQGKDM,2), NQ32, NSEG3, IOPTAU,  TCA0197
     *                LLAG, LLAGRS                                      TCA0197
      COMMON /LAGCM2/ TOLLAG
      COMMON /LAGCM3/ NPLAG
      DIMENSION IGS(10)                                                 TCA1196
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      LOGICAL LGS2(10)
      COMMON /LOGIC2/ LGS2
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     * CM2P
      COMMON /MORAB/  DAB, XKAB, AMAB, VDELTA
      COMMON /MORBC/  DBC, XKBC, AMBC
      COMMON /MORSTS/ DTS(NSADDM), XKTS(NSADDM), AMTS(NSADDM),
     * OMTS(NSADDM), OMIMG(NSADDM)
      COMMON /OPTION/ IOPT(20)
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
      COMMON /PSAGCM/ NPSAG
      COMMON /QUADKA/ PT(NQGKDM), WT(NQGKDM,2), NQ12, NSEG
      COMMON /QUADTH/ PT2(NQGKDM), WT2(NQGKDM,2), NQ22, NSEG2
      PARAMETER (NQWKDM=81)
      COMMON /QUADWK/ PTWK(NQWKDM), WTWK(NQWKDM), NQWK
      PARAMETER (NQCODM=81)
      COMMON /QUADCO/ PTCO(NQCODM), WTCO(NQCODM), NQCO
      LOGICAL LEQ
      COMMON /RATEQ/ RKQ(NTEMDM,2), RMVT(NTEMDM,2), RUS(NTEMDM,2),
     * LEQ(2)
      COMMON /RP2/    SLM, SLP
      COMMON /RP3/    NPRP
      COMMON /RP4/ STRUNL, STRUNR
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     * XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     * NSAD, NSADMX(2)
      COMMON /SIGCOM/ ISIGMA(2)                                         TCA1097
      COMMON /SINT/ SINT, PTS(NQGKDM), WTS(NQGKDM), NSS, NINTS
      COMMON /STATE/ TNP1, LSTATE, NSTATE
      COMMON /STATE2/ DGBND, LSBEND, NBND1, NBND2
      LOGICAL LSYM
      COMMON /SYM/ LSYM, NMID
      DIMENSION SCR(NQGKDM),EL(5,5),ENDPTS(2),NVAL(3),NEL(5),NDGEN(5,5) TCA1097
      EQUIVALENCE (ELFACT, EFACTF), (ELFACT(1,2), EFACTR)
      DATA GAM2T /1.3541179394D0/
      DATA AA /'Harmonic ', '  Morse  ', 'WKB-Morse'/
      DATA AB /' Harmonic ', 'Harm-Quart'/
C
C  constants
      PI = 4.0D0*ATAN(1.0D0)
      TPI = 2.0D0*PI
C  conversion from hartrees to eV
      EAU = 27.2113961D0                                                TCA1196
C  conversion from hartrees to kcal/mol
      CKCAL = 627.5095D0
C  conversion from hartrees to cm**-1
      CCM   = 219474.63067D0                                            TCA1196
C  beta = ckau/temp
      CKAU  = 315773.2D0                                                TCA1196
C  conversion from bohr radii to cm
      BOHR  = 0.529177249D-8                                            TCA1196
C  conversion from atomic unit of time to sec
      TAU   = 2.4188843d-17                                             TCA1196
C
C  all data input is in atomic units except where otherwise noted.
      WRITE (6, 600)
C  read title
      CALL TITLES (0, 5, 1)
C  write title
      CALL TITLES (1, 6, 1)
C  atomic labels
      READ (5, 502)  (AX(I), I=1,3)                                     TCA0197
      WRITE (6, 502) (AX(I), I=1,3)                                     TCA0197
C  atomic masses in amu
      READ (5, 500)  XMA, XMB, XMC
      WRITE (6, 500) XMA, XMB, XMC
C  The LGS and LGS2 variables are used internally.  The user will specify the
C  various options for the calculation by a series of integers IGS described
C  below.  The values of LGS and LGS2 variables are then set based on the
C  values of the IGS variables.
C
C  options:
C IGS( 1) - Allows the program to be stopped after certain stages.
C           1 - The program exits after computing reactant, product,
C               and saddle point properties, and TST rates.
C           2 - The program exits after the reaction coordinate and
C               the vibrationally adiabatic potential energy curve
C               are calculated.
C           3 - The program exits after the transmission coefficients
C               are calculated.
C           4 - Full calculation.
C IGS( 2) - Controls the type of calculation to be performed.
C           1 - Collinear calculation (1D) only.
C           2 - Three dimensional (3D) calculation only.
C           3 - Both 1D and 3D calculations are performed.
C           4 - Both 1D and 3D calculations are performed, but omit electronic
C               degeneracies and electronic excited states in the 1D
C               calculation.
C IGS( 3) - Selects the symmetry of the reaction path.
C           1 - The reaction path is symmetric about the saddle point.
C           2 - The reaction path is asymmetric (or is treated as asymmetric)
C               about the saddle point.
C IGS( 4) - Controls the reaction path calculation.
C           1 - The reaction path is calculated during program execution.
C           2 - The reaction path is read from LFN 1.
C           3 - The reaction path is calculated and written to LFN 1.
C IGS( 5) - Controls the treatment of the vibrational stretches.
C           1 - Harmonic method.
C           2 - Morse I method.
C           3 - WKB method.
C           4 - WKB method.  Results are written to LFN 14.
C           5 - WKB method.  Results are read from LFN 14.
C IGS( 6) - Controls treatment of the vibrational bends.
C           1 - Harmonic method.
C           2 - Harmonic-quartic method with energy levels computed using
C               the semiclassical method.
C           3 - Harmonic-quartic method with energy levels computed using
C               the semiclassical method if the harmonic force constant
C               is negative and using the perturbation-variation method
C               if the harmonic force constant is positive.
C           4 - Centrifugal oscillator method.
C IGS( 7) - Controls the tunneling methods used.
C           1 - Use MEPSAG and CD-SCSAG methods.
C           2 - Use MEPSAG and CD-SCSAG methods, and read LAG transmission
C               probabilities as needed from LFN 11 for 1D reactions, and
C               from LFN 12 for 3D reactions.
C           3 - Use MEPSAG, CD-SCSAG, LCG3, and LAG methods, but do not
C               perform a root search for the optimal alpha in the LAG
C               calulation.  The minimum in theta is taken from the three
C               points alpha = 0.0, 0.5, and 1.0.
C           4 - Use MEPSAG, CD-SCSAG, LCG3, and LAG methods.  The optimal
C               alpha is found in the LAG calculation.
C           5 - Use MEPSAG, CD-SCSAG, LCG3, and LAG methods, but do not
C               optimize alpha in the LAG calculation.  The minimum in theta
C               is taken from three points, alpha = 0.0, 0.5, and 1.0.  Write
C               the transmission probabilities from the large-curvature
C               tunneling methods (LCG3 and LAG) to LFN 11 for 1D reactions,
C               and to LFN 12 for 3D reactions.
C           6 - Use MEPSAG, CD-SCSAG, LCG3, and LAG methods, including
C               the optimization of alpha in the LAG calculation.  Write
C               the transmission probabilities from the large-curvature
C               tunneling methods (LCG3 and LAG) to LFN 11 for 1D reactions,
C               and to LFN 12 for 3D reactions.
C IGS( 8) - Controls printing of summary tables.
C           1 - Do not print any tables.
C           2 - Print a table of TST properties to LFN 21.
C           3 - Print a table of selected rate constants for LFN 22.
C           4 - Print both tables.
C IGS( 9) - Not used.
C IGS(10) - Not used.
C
C  LGS and LGS2 variables:
C     LGS(1) - if .true. print is on in datain
C              Printing is always on in v10.  User has no control.
C        (2) - if .true. execution is terminated after datain
C        (3) - if .true. morse correction to stretching vibration
C              is used.
C        (4) - if .true. quartic correction to bending vibration
C              used
C        (5) - if .true. execution is terminated after the
C              reaction coordinate is computed
C        (6) - if .true. execution is terminated after kapva
C        (7) - if .true. compute ground state eigenvalue by WKB
C        (8) - if .true. use semiclassical eigenvalues for bend
C              if .false. use semiclassical eigenvalues for bend if
C              fb<0, and perturbation-variation eigenvalues if fb>0
C        (9) - if .true. complete collinear calculation is omitted
C       (10) - if .true. complete 3d calculation is omitted
C    LGS2(1) - if .true. print bottleneck prop. of lfn 21
C        (2) - if .true. print brief summary on lfn 22
C        (3) - if .true., root search for optimum alpha not performed
C              in lag calculation
C        (4) - if .true., WKB energies, etc. write out to unit 14
C        (5) - if .true., WKB energies, etc. read in from unit 14
C        (6) - if .true., bend energy levels by centrifugal oscillator,
C              otherwise, treat as uncoupled degenerate oscillators.
C        (7) - if .true., RPH data written out to unit 1
C        (8) - if .true., LAG data written out to unit 11 (1D) 
C              and/or 12 (3D)
C        (9-10) - not used
C    IOPT(1) - if .ne. 0, state specific LAG calculation NSTATE (see
C              below) is initial state, and sum over final states from
C              IOPT(2) to IOPT(3)
C    IOPT(2) - if IOPT(1) .ne. 0, specifies final state
C    IOPT(3) - if IOPT(1) .ne. 0, specifies final state
C    IOPT(4) - if IOPT(1) .ne. 0, then the reactant asymptotic values
C              are not used in finding the maxima of VA and G.
C    IOPT(5) - if IOPT(1) .ne. 0, then the product asymptotic values
C              are not used in finding the maxima of VA and G.
C    IOPT(6) - controls printing of the free energy curves,
C              if IOPT(6) .eq. 0, then the free energy data are printed
C                 at all the temperatures for which the rate is calculated,
C              if IOPT(6) .eq. 1, then none of the free energy data are
C                 printed,
C              if IOPT(6) .eq. 2, then the free energy data are printed at
C                 the temperatures specified in the input.  When this option
C                 is specified, a line containing the number of temperatures
C                 at which the free energy data is to be printed is read
C                 immediately following the temperature data input.  Following
C                 this line is one or more lines of data containing the index
C                 of the temperature at which the free energy data is desired.
C    IOPT(7-20) - not used
      READ (5, 503) IGS, IOPT                                           TCA1196
      WRITE (6, 503) IGS, IOPT                                          TCA1196
C  Check validity of IGS input, and assign LGS and LGS2 variables,      TCA1196
C  and LLAG, LLAGRS, LRSTRT, and LSYM variables                         TCA1196
      LGS(1)=.TRUE.                                                     TCA1196
      IF ((IGS(1) .LT. 1) .OR. (IGS(1) .GT. 4)) THEN                    TCA1196
        WRITE(6,6001) 1                                                 TCA1196
        STOP                                                            TCA1196
      ELSE IF (IGS(1) .EQ. 1) THEN                                      TCA1196
        LGS(2)=.TRUE.                                                   TCA1196
        LGS(5)=.FALSE.                                                  TCA1196
        LGS(6)=.FALSE.                                                  TCA1196
      ELSE IF (IGS(1) .EQ. 2) THEN                                      TCA1196
        LGS(2)=.FALSE.                                                  TCA1196
        LGS(5)=.TRUE.                                                   TCA1196
        LGS(6)=.FALSE.                                                  TCA1196
      ELSE IF (IGS(1) .EQ. 3) THEN                                      TCA1196
        LGS(2)=.FALSE.                                                  TCA1196
        LGS(5)=.FALSE.                                                  TCA1196
        LGS(6)=.TRUE.                                                   TCA1196
      ELSE IF (IGS(1) .EQ. 4) THEN                                      TCA1196
        LGS(2)=.FALSE.                                                  TCA1196
        LGS(5)=.FALSE.                                                  TCA1196
        LGS(6)=.FALSE.                                                  TCA1196
      END IF                                                            TCA1196
      IF ((IGS(2) .LT. 1) .OR. (IGS(2) .GT. 4)) THEN                    TCA0897
        WRITE(6,6001) 2                                                 TCA1196
        STOP                                                            TCA1196
      ELSE IF (IGS(2) .EQ. 1) THEN                                      TCA1196
        LGS(9)=.FALSE.                                                  TCA1196
        LGS(10)=.TRUE.                                                  TCA1196
      ELSE IF (IGS(2) .EQ. 2) THEN                                      TCA1196
        LGS(9)=.TRUE.                                                   TCA1196
        LGS(10)=.FALSE.                                                 TCA1196
      ELSE IF ((IGS(2) .EQ. 3) .OR. (IGS(2) .EQ. 4)) THEN               TCA0897
        LGS(9)=.FALSE.                                                  TCA1196
        LGS(10)=.FALSE.                                                 TCA1196
      END IF                                                            TCA1196
      LGS(9) = .NOT.LGS(9)
      LGS(10) = .NOT.LGS(10)
      IF ((IGS(3) .LT. 1) .OR. (IGS(3) .GT. 2)) THEN                    TCA1196
        WRITE(6,6001) 3                                                 TCA1196
        STOP                                                            TCA1196
      ELSE IF (IGS(3) .EQ. 1) THEN                                      TCA1196
        LSYM=.TRUE.                                                     TCA1196
      ELSE IF (IGS(3) .EQ. 2) THEN                                      TCA1196
        LSYM=.FALSE.                                                    TCA1196
      END IF                                                            TCA1196
      IF ((IGS(4) .LT. 1) .OR. (IGS(4) .GT. 3)) THEN                    TCA1196
        WRITE(6,6001) 4                                                 TCA1196
        STOP                                                            TCA1196
      ELSE IF (IGS(4) .EQ. 1) THEN                                      TCA1196
        LGS2(7)=.FALSE.                                                 TCA1196
        LRSTRT=.FALSE.                                                  TCA1196
      ELSE IF (IGS(4) .EQ. 2) THEN                                      TCA1196
        LGS2(7)=.FALSE.                                                 TCA1196
        LRSTRT=.TRUE.                                                   TCA1196
      ELSE IF (IGS(4) .EQ. 3) THEN                                      TCA1196
        LGS2(7)=.TRUE.                                                  TCA1196
        LRSTRT=.FALSE.                                                  TCA1196
      END IF                                                            TCA1196
      IF ((IGS(5) .LT. 1) .OR. (IGS(5) .GT. 5)) THEN                    TCA1196
        WRITE(6,6001) 5                                                 TCA1196
        STOP                                                            TCA1196
      ELSE IF (IGS(5) .EQ. 1) THEN                                      TCA1196
        LGS(3)=.FALSE.                                                  TCA1196
        LGS(7)=.FALSE.                                                  TCA1196
        LGS2(4)=.FALSE.                                                 TCA1196
        LGS2(5)=.FALSE.                                                 TCA1196
      ELSE IF (IGS(5) .EQ. 2) THEN                                      TCA1196
        LGS(3)=.TRUE.                                                   TCA1196
        LGS(7)=.FALSE.                                                  TCA1196
        LGS2(4)=.FALSE.                                                 TCA1196
        LGS2(5)=.FALSE.                                                 TCA1196
      ELSE IF (IGS(5) .EQ. 3) THEN                                      TCA1196
        LGS(3)=.TRUE.                                                   TCA1196
        LGS(7)=.TRUE.                                                   TCA1196
        LGS2(4)=.FALSE.                                                 TCA1196
        LGS2(5)=.FALSE.                                                 TCA1196
      ELSE IF (IGS(5) .EQ. 4) THEN                                      TCA1196
        LGS(3)=.TRUE.                                                   TCA1196
        LGS(7)=.TRUE.                                                   TCA1196
        LGS2(4)=.TRUE.                                                  TCA1196
        LGS2(5)=.FALSE.                                                 TCA1196
      ELSE IF (IGS(5) .EQ. 5) THEN                                      TCA1196
        LGS(3)=.TRUE.                                                   TCA1196
        LGS(7)=.TRUE.                                                   TCA1196
        LGS2(4)=.FALSE.                                                 TCA1196
        LGS2(5)=.TRUE.                                                  TCA1196
      END IF                                                            TCA1196
      IF ((IGS(6) .LT. 1) .OR. (IGS(6) .GT. 4)) THEN                    TCA1196
        WRITE(6,6001) 6                                                 TCA1196
        STOP                                                            TCA1196
      ELSE IF (IGS(6) .EQ. 1) THEN                                      TCA1196
        LGS(4)=.FALSE.                                                  TCA1196
        LGS(8)=.FALSE.                                                  TCA1196
        LGS2(6)=.FALSE.                                                 TCA1196
      ELSE IF (IGS(6) .EQ. 2) THEN                                      TCA1196
        LGS(4)=.TRUE.                                                   TCA1196
        LGS(8)=.TRUE.                                                   TCA1196
        LGS2(6)=.FALSE.                                                 TCA1196
      ELSE IF (IGS(6) .EQ. 3) THEN                                      TCA1196
        LGS(4)=.TRUE.                                                   TCA1196
        LGS(8)=.FALSE.                                                  TCA1196
        LGS2(6)=.FALSE.                                                 TCA1196
      ELSE IF (IGS(6) .EQ. 4) THEN                                      TCA1196
        LGS(4)=.TRUE.                                                   TCA0497
        LGS(8)=.FALSE.                                                  TCA1196
        LGS2(6)=.TRUE.                                                  TCA1196
      END IF                                                            TCA1196
      IF ((IGS(7) .LT. 1) .OR. (IGS(7) .GT. 6)) THEN                    TCA1196
        WRITE(6,6001) 7                                                 TCA1196
        STOP                                                            TCA1196
      ELSE IF (IGS(7) .EQ. 1) THEN                                      TCA1196
        LGS2(3)=.FALSE.                                                 TCA1196
        LGS2(8)=.FALSE.                                                 TCA1196
        LLAG=.FALSE.                                                    TCA1196
        LLAGRS=.FALSE.                                                  TCA1196
      ELSE IF (IGS(7) .EQ. 2) THEN                                      TCA1196
        LGS2(3)=.FALSE.                                                 TCA1196
        LGS2(8)=.FALSE.                                                 TCA1196
        LLAG=.TRUE.                                                     TCA1196
        LLAGRS=.TRUE.                                                   TCA1196
      ELSE IF (IGS(7) .EQ. 3) THEN                                      TCA1196
        LGS2(3)=.TRUE.                                                  TCA1196
        LGS2(8)=.FALSE.                                                 TCA1196
        LLAG=.TRUE.                                                     TCA1196
        LLAGRS=.FALSE.                                                  TCA1196
      ELSE IF (IGS(7) .EQ. 4) THEN                                      TCA1196
        LGS2(3)=.FALSE.                                                 TCA1196
        LGS2(8)=.FALSE.                                                 TCA1196
        LLAG=.TRUE.                                                     TCA1196
        LLAGRS=.FALSE.                                                  TCA1196
      ELSE IF (IGS(7) .EQ. 5) THEN                                      TCA1196
        LGS2(3)=.TRUE.                                                  TCA1196
        LGS2(8)=.TRUE.                                                  TCA1196
        LLAG=.TRUE.                                                     TCA1196
        LLAGRS=.FALSE.                                                  TCA1196
      ELSE IF (IGS(7) .EQ. 6) THEN                                      TCA1196
        LGS2(3)=.FALSE.                                                 TCA1196
        LGS2(8)=.TRUE.                                                  TCA1196
        LLAG=.TRUE.                                                     TCA1196
        LLAGRS=.FALSE.                                                  TCA1196
      END IF                                                            TCA1196
      IF ((IGS(8) .LT. 1) .OR. (IGS(8) .GT. 4)) THEN                    TCA1196
        WRITE(6,6001) 8                                                 TCA1196
        STOP                                                            TCA1196
      ELSE IF (IGS(8) .EQ. 1) THEN                                      TCA1196
        LGS2(1)=.FALSE.                                                 TCA1196
        LGS2(2)=.FALSE.                                                 TCA1196
      ELSE IF (IGS(8) .EQ. 2) THEN                                      TCA1196
        LGS2(1)=.TRUE.                                                  TCA1196
        LGS2(2)=.FALSE.                                                 TCA1196
      ELSE IF (IGS(8) .EQ. 3) THEN                                      TCA1196
        LGS2(1)=.FALSE.                                                 TCA1196
        LGS2(2)=.TRUE.                                                  TCA1196
      ELSE IF (IGS(8) .EQ. 4) THEN                                      TCA1196
        LGS2(1)=.TRUE.                                                  TCA1196
        LGS2(2)=.TRUE.                                                  TCA1196
      END IF                                                            TCA1196
C
C  number of quadrature points for kappa and theta
C  lag options added to next card
C  quadrature info for s integration in LAG2 add also
      READ (5, 510) NSEG, NQ12, NSEG2, NQ22, NSEG3, NQ32, IOPTAU,       TCA1196
     * SINT, NSS, NINTS, TOLLAG                                         TCA0997
      WRITE (6, 510) NSEG, NQ12, NSEG2, NQ22, NSEG3, NQ32, IOPTAU,      TCA1196
     * SINT, NSS, NINTS, TOLLAG                                         TCA0997
C  order of interpolation used for reaction coordinate values,
C  the order of Legendre quadrature for WKB phase integrals and
C  centrifugal oscillator phase integrals, and the maximum primary 
C  quantum number in partition functions of cent. osc.  (Note that
C  the h.o. energy levels are (v+1)hbar*omega = (2*n+K+1)hbar*omega
C  and NLEVEL is maximum in v not n)
      READ (5, 504)  NINT, NQWK, NQCO, NLEVEL
      WRITE (6, 504) NINT, NQWK, NQCO, NLEVEL
C  approximate asymptotic geometries
      READ (5, 500)  R1ASY, R2ASY
      WRITE (6, 500) R1ASY, R2ASY
C  number of saddle points
      READ (5, 504)  NSAD
      WRITE (6, 504) NSAD
      LSAD = NSAD.GT.0
      N = NSAD
      IF (.NOT.LSAD) N = 1
C  approximate saddle point geometries
      READ (5, 501)  (R1SP(I), R2SP(I), I=1,N)
      WRITE (6, 501) (R1SP(I), R2SP(I), I=1,N)
C  symmetry and restart options, parameters for extrapolating the
C  curvature and truncating the reaction path.  SM2,SP2 parameters used
C  to extrapolate curvature; they are values of reaction coordinate at
C  which curvature is set to half the value at SM,S, respectively.  If
C  SM2(SP2) is less than SM(S) then no extrpolation is done.  STRUNL,
C  STRUNR are values of reaction coordinate at which the reaction path
C  is truncated on the left and right, respectively.  If both are zero
C  no truncation is performed.
      READ (5, 505)   SM2, SP2, STRUNL, STRUNR
      WRITE (6, 505)  SM2, SP2, STRUNL, STRUNR
      IF (.NOT.LRSTRT) THEN
C  parameters for reaction coordinate calculation
C  bounds on reaction coordinate calculation
         READ (5, 500)  SLM, SLP
         WRITE (6, 500) SLM, SLP
C  DEL - step size, DELSV - interval along s at which potential
C  parameter are stored, ACASY - if acos(diff.  between the gradient and
C  the asymptotic value of the grad) is .gt.  ACASY then it is assumed
C  in the asymptotic region, EPSASY - if relative diff.  between the
C  potential parameters at S and the asymptotic values are less than
C  EPSASY reaction coor.  calc.  is ended
         READ (5, 500)  DEL, DELSV, ACASY, EPSASY
         WRITE (6, 500) DEL, DELSV, ACASY, EPSASY
C  parameters for the calculation of the curvature.  DELCUR - step size
C  central difference derivative, DELCSP - special value of DELCUR near
C  the saddle point, SM,SP - bounds on values of s at which curvature is
C  calculated
         READ (5, 500)  DELCUR, DELCSP, SM, SP
         WRITE (6, 500) DELCUR, DELCSP, SM, SP
C  option for fitting the bending potential
C  IBOPT=3, fit using two angles, THETA(1) and THETA(2)
C  IBOPT.LT.3, FB fit to second derivative at THETA=0
C  if FB.GT.0 AB fit to fourth derivative at THETA=0
C  if FB.LT.0, IBOPT=1, fit FB,AB to give position and depth of well
C              IBOPT=2, fit FB,AB to give well depth and tp for E=V0
C  state selected option - LSTATE and vibrational state
C  limits in s for searching for local maximum in adiabatic barrier
C  bending state selection
C  state selected option for bend - LSBEND and bending quantum numbers
C     If LGS2(6)=.F., bends are treated as uncoupled oscillators and
C     NBND1,NBND2 are the quantum numbers of the independent modes.
C     If LGS2(6)=.T., bends are treated as coupled by the centrifugal
C     oscillator method and NBND1=n, NBND2=K where the harmonic energy
C     levels are (2n+K+1)*hbar*omega
      END IF
      READ (5,508) IBOPT, THETA, LSTATE, NSTATE, SAD1, SAD2, LSBEND,
     * NBND1, NBND2
      WRITE (6,508) IBOPT, THETA, LSTATE, NSTATE, SAD1, SAD2, LSBEND,
     * NBND1, NBND2
      IF (LSTATE .EQ. 0) NSTATE = 0
      TNP1 = 2.D0*DBLE(NSTATE) + 1.D0                                   GCL1092
      IF (LSBEND .EQ. 0) NBND1 = 0
      IF (LSBEND .EQ. 0) NBND2 = 0
      DGBND = 1.0D0
      IF (LGS2(6)) THEN
         IF(NBND2.EQ.0) DGBND = 2.0D0
      ELSE
         IF (NBND1.NE.NBND2) DGBND = 2.0D0
      END IF
C  number of temperatures
C  NPGTST - print interval for free energy curves. 0 - print all
C  NPRP   - print interval for reaction path summary. 0 - print all
C  NPAD   - print interval for adiabatic potential curves. 0 - print all
C  NPSAG  - print option for sag calculations. 0 - print all, all else
C           suppress printing probabilities
C  NPLAG  - print option for lag calculations. 0 - print all, all else
C           suppress printing probabilities
      READ (5, 504)  NTMAX, NPGTST, NPRP, NPAD, NPSAG, NPLAG
      WRITE (6, 504) NTMAX, NPGTST, NPRP, NPAD, NPSAG, NPLAG
C  check that NTMAX is in [1,100]
      IF ((NTMAX .LT. 1) .OR. (NTMAX .GT. 100))                         TCA0197
     *  STOP 'NTMAX is not in the range [1,100]'                        TCA0197
C  read temperatures
      READ (5, 530)  (TEMP(IT), IT=1, NTMAX)                            TCA0197
      WRITE (6, 530) (TEMP(IT), IT=1, NTMAX)                            TCA0197
C  option to read in other coll. rates and compare to calc'd rates
      READ (5, 509)  (LEQ(I), I=1,2)
      WRITE (6, 509) (LEQ(I), I=1,2)
      DO 15 IC3D = 1,2
         IF (.NOT.LEQ(IC3D)) GO TO 15
C  exact rate constants at temperatures read in above
         READ (5, 507)  (RKQ(I, IC3D), I=1,NTMAX)
         WRITE (6, 507) (RKQ(I, IC3D), I=1,NTMAX)
C  muvt rate constants
         READ (5, 507)  (RMVT(I, IC3D), I=1,NTMAX)
         WRITE (6, 507) (RMVT(I, IC3D), I=1,NTMAX)
C  unified statistical (microcan.) rate constants
         READ (5, 507)  (RUS(I, IC3D), I=1,NTMAX)
         WRITE (6, 507) (RUS(I, IC3D), I=1,NTMAX)
15    CONTINUE
C  prepare for IOPT(6) options...
C  first initialize the PTEMP array
      DO 10 I=1,NTEMDM                                                  TCA0996
        IF (IOPT(6) .EQ. 0) THEN                                        TCA0996
          PTEMP(I)=.TRUE.                                               TCA0996
        ELSE                                                            TCA0996
          PTEMP(I)=.FALSE.                                              TCA0996
        END IF                                                          TCA0996
10    CONTINUE                                                          TCA0996
C  if IOPT(6) .eq. 2, then we need to read in the indices of the temperatures
C  at which the free energy data is to be printed.
C  first read a single line containing the number of indices to read
      IF (IOPT(6) .EQ. 2) THEN                                          TCA0996
        READ (5,520) NIDX                                               TCA0996
        WRITE (6,520) NIDX                                              TCA0996
C  now read in the temperature indices
        READ (5,521) (ITEMP(IT), IT=1,NIDX)                             TCA0996
        WRITE (6,521) (ITEMP(IT), IT=1,NIDX)                            TCA0996
C  and set up the PTEMP array so we know when to print
        DO 11 I=1,NIDX                                                  TCA0996
          IF (ITEMP(I) .LT. 1 .OR. ITEMP(I) .GT. NTMAX) THEN            TCA0996
            STOP 'Error in datain - invalid temperature index.'         TCA0996
          END IF                                                        TCA0996
          PTEMP(ITEMP(I))=.TRUE.                                        TCA0996
11      CONTINUE                                                        TCA0996
      END IF                                                            TCA0996
C  number of activation energies to be computed
      READ (5, 504)  IACT
      WRITE (6, 504) IACT
C  indices of temperatures used in act. energy calc.
      READ (5, 504)  (NT1(IT), NT2(IT), IT=1,IACT)
      WRITE (6, 504) (NT1(IT), NT2(IT), IT=1,IACT)
C  electronic partition function parameters
C  number of electronic energy levels for A,BC,ABC,AB,C
      READ (5, 504)  (NEL(I), I=1,5)
      WRITE (6, 504) (NEL(I), I=1,5)
      DO 20 I = 1,5
         EL(1,I) = 0.D0                                                 GCL1092
         NDGEN(1,I) = 1
         IF (NEL(I) .EQ. 0) GO TO 20
         N = NEL(I)
C  degeneracies of the energy level
         READ (5, 504)  (NDGEN(J, I), J=1,N)
         WRITE (6, 504) (NDGEN(J, I), J=1,N)
C  energy levels in kcal/mol
         READ (5, 506)  (EL(J, I), J=1,N)
         WRITE (6, 506) (EL(J, I), J=1,N)
         NEL(I) = MAX(1,NEL(I))
20    CONTINUE
C  special energies for tunneling calculations
C  total energies in kcal/mol
      NESPEC = 0
      READ (5, 504, END=40)  NESPEC
      WRITE (6, 504) NESPEC
      IF (NESPEC .GT. 0) THEN
         READ (5, 500)  (ESPEC(I), I=1,NESPEC)
         WRITE (6, 500) (ESPEC(I), I=1,NESPEC)
         DO 30 I = 1,NESPEC
            ESPEC(I) = ESPEC(I)/CKCAL
30       CONTINUE
      END IF
40    CONTINUE
C
C   Check that all restart files exist
      CALL CHKRST                                                       GCL1096
C
      CALL PREPEF                                                       GCL0992
C
C  OPTIONS
      WRITE (6, 606)
      IF (LGS(9) .AND. LGS(10)) THEN                                    GCL0394
          WRITE (6, 609)                                                GCL0394
      ELSE IF (LGS(9)) THEN                                             GCL0394
          WRITE (6, 608)                                                GCL0394
      ELSE                                                              GCL0394
          WRITE (6, 610)                                                GCL0394
      ENDIF                                                             GCL0394
      IF (LGS(3)) WRITE (6, 612)
      IF (.NOT.LGS(3)) WRITE (6, 614)
      IF (LGS(7)) THEN
         CALL GAUSSQ(1, NQWK, 0.D0, 0.D0, 0, ENDPTS, SCR, PTWK, WTWK)
         WRITE (6, 616) NQWK
         IF (LGS2(5)) THEN
            WRITE (6, 617)
            CALL WKBRD
         END IF
      END IF
      IF (IBOPT .LT. 1) IBOPT=1
      IF (IBOPT .GT. 3) IBOPT=3
      IF (IBOPT .EQ. 1) WRITE (6, 620)
      IF (IBOPT .EQ. 2) WRITE (6, 622)
      IF (IBOPT .EQ. 3) WRITE (6, 624) THETA
      IF (LGS2(6)) THEN
         CALL GAUSSQ(1, NQCO, 0.D0, 0.D0, 0, ENDPTS, SCR, PTCO, WTCO)
         WRITE (6, 619) NQCO, NLEVEL
      ELSE IF (LGS(4)) THEN
         WRITE (6, 618)
         IF (LGS(8)) WRITE (6, 626)
         IF (.NOT.LGS(8)) WRITE (6, 628)
      ELSE
         WRITE (6, 630)
      END IF
      IF (LSTATE .NE. 0) WRITE (6, 632) NSTATE
      IF (LSBEND .NE. 0) WRITE (6, 634) NBND1,NBND2
      IF (IOPT(1) .EQ. 0) THEN
         IOPT(2) = NSTATE
         IOPT(3) = NSTATE
      ELSE
         IF (IOPT(2) .LT. NSTATE) THEN
            WRITE (6, 682) IOPT(2), NSTATE
            STOP 'DATAIN 2'
         END IF
         IF (IOPT(3) .LT. IOPT(2)) THEN
            WRITE (6, 684) IOPT(2), IOPT(3)
            STOP 'DATAIN 3'
         END IF
      END IF
      THETA(1) = PI*THETA(1)/180.D0                                     GCL1092
      THETA(2) = PI*THETA(2)/180.D0                                     GCL1092
      A12 = THETA(1)*THETA(1)
      A11 = 0.5D0*A12
      A12 = A12*A12/24.D0                                               GCL1092
      A22 = THETA(2)*THETA(2)
      A21 = 0.5D0*A22
      A22 = A22*A22/24.D0                                               GCL1092
      THETA(1) = PI - THETA(1)
      THETA(2) = PI - THETA(2)
      IF (LLAG) THEN
         WRITE (6, 636)
         IF (IOPT(1) .NE. 0) THEN
            WRITE (6, 638) IOPT(2), IOPT(3)
         END IF
         IF (IOPTAU .EQ. 0) WRITE (6, 640)
         IF (IOPTAU .NE. 0) WRITE (6, 642)
         WRITE (6, 646)
         WRITE (6, 648) SINT, NSS, NINTS
         IF (LGS2(3)) THEN
            WRITE (6, 650)
         ELSE
            WRITE (6, 652) TOLLAG
         END IF
      END IF
C
C  ISIGMA(1) and (2) are symmetry numbers for the rotational
C     partition functions for reactants and products,
C     respectively.
C  ISIGMA(1) and (2) are also the reaction coordinate degeneracies
C     for the forward and reverse reactions, respectively.
      ISIGMA(1) = 1
      ISIGMA(2) = 1
      IF (XMA .EQ. XMB) ISIGMA(2) = 2
      IF (XMB .EQ. XMC) ISIGMA(1) = 2
C  masses
      IF (LGS(1)) WRITE (6, 654) AX(1), XMA, AX(2), XMB, AX(3), XMC
      XM = XMB + XMC
C     center of mass of BC pair
      CM1 = XMC/XM
      XMU = XMA*XM/(XMA+XM)
      XM = XMB*XMC/XM
C     scaling parameter
      CM2 = SQRT(XM/XMU)
      XMP = XMA + XMB
      CM1P = XMA/XMP
      XMUP = XMC*XMP/(XMC+XMP)
      XMP = XMA*XMB/XMP
      CM2P = SQRT(XMP/XMUP)
C     skew angle
      T = CM2/CM1
      T = ATAN(T)*180.D0/PI
      IF (LGS(1)) WRITE (6, 656) AX(2), AX(3), XM, AX(1), AX(2), XMP,
     * AX(1), AX(2), AX(3), XMU, CM1, CM2, T
C
C  DBC - value of the potential for all three atoms infinitly separated
C     search for DBC
      R1 = 15.D0                                                        GCL1092
      R2 = 15.D0                                                        GCL1092
      V = 0.D0                                                          GCL1092
50    CONTINUE
      V0 = V
      R1 = R1+5.D0                                                      GCL1092
      R2 = R2+5.D0                                                      GCL1092
      R3 = R1+R2
      CALL PEF (R1, R2, R3, V, D1, D2, D3, 0)                           GCL0893
      IF (ABS((V-V0)/V) .GT. 1.D-10) GO TO 50
      DBC = V
C
C  zero of energy taken to be the bottom of the reactant diatomic well.
C     search for asymptotic equilibrium geometries
      CALL GEOM (LGS(1), AX(2), AX(3), 0, DBC, R2ASY, V0, DBC, OMBC,
     * XKBC, AMBC)
      IF (V0 .GE. 1.D-6) THEN
         WRITE (6, 6002) V0
         STOP 'DATAIN 5'
      END IF
      CALL GEOM (LGS(1), AX(1), AX(2), 1, DBC, R1ASY, VDELTA, DAB,
     * OMAB, XKAB, AMAB)
      IF (LGS(1)) THEN
         T1 = DBC*EAU
         T2 = DBC*CKCAL
         T3 = VDELTA*EAU
         T4 = VDELTA*CKCAL
         WRITE (6, 658) DBC, T1, T2, VDELTA, T3, T4
      END IF
C
C  saddle point calculation
      IF (.NOT.LSAD) THEN
C     no saddle point, saddle point parameters use asymptotic values
         WRITE (6, 660)
         IF (VDELTA .GE. 0.D0) THEN
            WRITE (6, 6003)
            STOP 'DATAIN 6'
         END IF
         R2SP(1) = R2ASY
         XSP(1) = R1SP(1) + CM1*R2SP(1)
         YSP(1) = CM2*R2SP(1)
         DX = 1.D0                                                      GCL1092
         DY = 0.D0                                                      GCL1092
         CALL VMIN (XSP(1), YSP(1), DX, DY, 0.1D0, VSP(1), .FALSE.,
     *    IERR)
         IF (IERR .NE. 0) THEN
            WRITE (6, 6004)
            STOP 'DATAIN 7'
         END IF
         SVECT(1,1) = -1.D0                                             GCL1092
         SVECT(2,1) = 0.D0                                              GCL1092
         UVECT(1,1) = 0.D0                                              GCL1092
         UVECT(2,1) = 1.D0                                              GCL1092
         CALL D2VDU2 (XSP(1), YSP(1), 0.D0, 1.D0, D2V)
         DTS(1) = DBC - VSP(1)
         OMTS(1) = SQRT(D2V/XMU)
         XKTS(1) = 4.D0*DTS(1)/OMTS(1)
         AMTS(1) = 2.D0*SQRT(2.D0*XM*DTS(1))/XKTS(1)
         OMIMG(1) = 0.D0                                                GCL1092
         CALL BEND (XSP(1), YSP(1), FBTS(1), QFBTS(1), GBTS(1), XMOMTS
     *    (1))
         NSADMX(1) = 1
         NSADMX(2) = 1
      ELSE
C     root search for saddle point configuration
         CALL SADDLE
      END IF
C
C  electronic factors
      CALL ELPART (AX, NEL, EL, NDGEN)
      IF (IGS(2) .EQ. 4) WRITE(6,661)                                   TCA0897
C
C  partition functions
      IF (LGS(3)) THEN
         IF (LGS(7)) THEN
            IMH = 3
         ELSE
            IMH = 2
         END IF
      ELSE
         IMH = 1
      END IF
      IF (LGS(1)) THEN
         WRITE (6, 662) AA(IMH)
         IF (LSTATE .NE. 0) WRITE (6, 664) NSTATE
         WRITE (6, 666) AA(IMH)
      END IF
      XIP = XMP*R1ASY*R1ASY
      XI = XM*R2ASY*R2ASY
      DO 60 IT = 1,NTMAX
         BET = BETA(IT)
         RT = CKCAL/BET
         T = SQRT(XMU/(TPI*BET))/BOHR
         TP = SQRT(XMUP/(TPI*BET))/BOHR
         QTRCOL = LOG(T)
         QTR = 3.D0*QTRCOL
         TP = LOG(TP)
         QTRP = 3.D0*TP
C  all partition function subroutines return log(pf)
         CALL PFCNST (BET, -1, -1000.D0, DBC, XKBC, QSTRM, QSC, .TRUE.,
     *    LGS(7))
         RATIO(IT,1) = QSC
         CALL PFCNST (BET, -2, 1000.D0, DAB, XKAB, QSTRP, QSCP, .TRUE.,
     *    LGS(7))
         QSTR = QSTRM
         CALL PFCNR (BET, XI, ISIGMA(1), QROTCL, QROTQ)                 GCL0795
         QRC = EXP(QROTQ - QROTCL)
         RATIO(IT,2) = QRC
         CALL PFCNR (BET, XIP, ISIGMA(2), QRTCLP, QROTQP)               GCL0795
         T = QSTR - EFACTF(IT)
         CPHIC(IT) = QTRCOL + T
         IF (IGS(2) .EQ. 4) CPHIC(IT) = CPHIC(IT) + EFACTF(IT)          TCA0897
         CPHI(IT) = QTR + QROTQ + T
         T = QSTRP - EFACTR(IT)
         EQUIL(IT,1) = TP + T - CPHIC(IT)
         IF (IGS(2) .EQ. 4) EQUIL(IT,1) = EQUIL(IT,1) + EFACTR(IT)      TCA0897
         EQUIL(IT,2) = QTRP + T + QROTQP - CPHI(IT)
         T = BET*VDELTA
         EQUIL(IT,1) = EQUIL(IT,1) - T
         EQUIL(IT,2) = EQUIL(IT,2) - T
         CNST(IT) = RT*CPHI(IT)
         CNSTC(IT) = RT*CPHIC(IT)
         IF (LGS(1)) WRITE (6, 668) TEMP(IT), BET, QTRCOL, QTR, QROTQ,
     *    ISIGMA(1), QRC, QSTR, QSC, CPHI(IT), EQUIL(IT,1),
     *    EQUIL(IT,2)
60    CONTINUE
C
      DO 80 IC3D = 1,2
         IF (.NOT.LGS(IC3D+8)) GO TO 80
         NMX = NSADMX(IC3D)
         WSTR = OMIMG(NMX)
         DO 70 IT = 1,NTMAX
            T = BETA(IT)*WSTR
            KAPW(IT,IC3D) = 1.D0 + T*T/24.D0                            GCL1092
70       CONTINUE
80    CONTINUE
      IF (LGS(1) .AND. LSAD) THEN
         IF (LGS(4)) THEN
            IHQ = 2
         ELSE
            IHQ = 1
         END IF
         WRITE (6, 670) AA(IMH), AB(IHQ)
         IF (LSTATE .NE. 0) WRITE (6, 672) NSTATE
         DO 100 IC3D = 1,2
            IF (.NOT.LGS(IC3D+8)) GO TO 100
            IF (IC3D .EQ. 1) THEN
               QB = 0.D0                                                GCL1092
               QROTQ = 0.D0                                             GCL1092
               CB = 0.D0                                                GCL1092
               CROT = 0.D0                                              GCL1092
               WRITE (6, 673)
            ELSE
               WRITE (6, 674)
            END IF
            NMX = NSADMX(IC3D)
            WRITE (6, 676) AA(IMH), AB(IHQ)
            DO 90 IT = 1,NTMAX
               BET = BETA(IT)
               CCC = EXP(-BET*VSP(NMX))/(TPI*BET*TAU)
               IS = -NMX-2
               CALL PFCNST (BET, IS, 0.D0, DTS(NMX), XKTS(NMX), QSTR,
     *          CSTR, .TRUE., LGS(7))
               IF (IC3D .EQ. 2) THEN
                  CALL PFCNB(BET, FBTS(NMX), QFBTS(NMX), GBTS(NMX), QB,
     *             CB, DELPHI, .TRUE.)
                  CALL PFCNR(BET, XMOMTS(NMX), ISIGMA(1), QROTCL,       TCA1097
     *             QROTQ)                                               GCL0795
                  CROT = EXP(QROTQ - QROTCL)
               END IF
               WC = KAPW(IT,IC3D)
               IF (IC3D .EQ. 1) THEN
                  RK = CCC*EXP(QSTR-CPHIC(IT))
               ELSE
                  RK = CCC*EXP(QSTR+QB+QROTQ-CPHI(IT))*DBLE(ISIGMA(1))  TCA1097
               END IF
               WRITE (6, 678) TEMP(IT), QSTR, CSTR, QB, CB, QROTQ,
     *          ISIGMA(1), CROT, WC, RK, RK*WC
90          CONTINUE
100      CONTINUE
      END IF
      RETURN
500   FORMAT (4F20.10)
501   FORMAT (2F20.10)
502   FORMAT (5X, 3(3X,A2))                                             TCA0197
503   FORMAT (10I2, 2X, 20I2)                                           TCA1196
504   FORMAT (10I5)
505   FORMAT (2F20.10, 2F15.8)                                          TCA1196
506   FORMAT (8F10.4)
507   FORMAT (4G20.10)
508   FORMAT (I5, 2F10.4, 1X, I1, I5, 2F20.10, 1X, I1, 2I3)
509   FORMAT (1X, L1, 5X, L1)                                           TCA0997
510   FORMAT (7I5, F10.5, 2I5, E10.2)                                   TCA0997
520   FORMAT (2X,I3)                                                    TCA0996
521   FORMAT (10(I5))                                                   TCA0996
530   FORMAT (10(1X,F7.1))                                              TCA0197
600   FORMAT (/,1X,T5,'Transition state rate constants, reactants ',
     *   'are A + BC',//,1X,T5,'Input data',/)
606   FORMAT (/,1X,34('*'),2X,'Options',1X,34('*'))
608   FORMAT (/,1X,T2,'The collinear calculation will be done.')
609   FORMAT (/,1X,T2,'The collinear and 3D calculations will be done.')
610   FORMAT (/,1X,T2,'The 3D calculation will be done.')
612   FORMAT (/,1X,T2,'Vibrational energy levels:',/,1X,T5,
     * 'Morse approximation I will be used for the stretch.')
614   FORMAT (1X,T2,'Vibrational energy levels:',/,1X,T5,
     * 'Harmonic approximation will be used for the stretch.')
616   FORMAT (1X,T5,'The ground vibrational eigenvalue for the',
     * ' stretch is determined by',/,1X,T5,'the WKB approximation.',
     *        /,1X,T5,'For state-selected rates the NSTR level',
     * ' is also determined by',/,1X,T5,'the WKB approximation.',
     *        /,1X,T5,'Phase integrals for WKB eigenvalues are ',
     * 'found by a ', I4, ' point',/,1X,T5,'Legendre quadrature.')
617   FORMAT (1X,T5,'The WKB energy levels are read in from unit 14.')
618   FORMAT (/,1X,T5,'Harmonic-quartic approximation is used for ',
     *                'the bend.')
619   FORMAT (/,1X,T5, 'The bend energy levels are computed by ',
     *                 'the centrifugal oscillator.',
     *        /,1X,T5, 'The phase integrals are computed by ',
     *             I4, ' point Legendre quadrature.',
     *        /,1X,T5, 'The bend partition functions are ',
     *                 'obtained by summing over the energy levels ',
     *        /,1X,T5, 'up to v2 = ', I4, 
     *                 ' where K = v2, v2-2,...,1 or 0 for each v2.')
620   FORMAT (1X,T5,'IBOPT=1, the harmonic force constant FB is ',
     * 'fitted using the second',/,1X,T5,'derivative w.r.t. the ',
     * 'angle at the equilibrium collinear geometry.',
     * /,1X,T5,'If FB is positive the quartic force constant is ',
     * 'fitted using the',/,1X,T5,'fourth derivative w.r.t. the ',
     * 'angle.  If FB is negative the harmonic',/,1X,T5,'and quartic ',
     * 'force constants are fitted to reproduce the well depth',
     * /,1X,T5,'and the turning point for an energy equal to the ',
     * 'local maximum of',/,1X,T5,'the double well.')
622   FORMAT (1X,T5,'IBOPT=2, the harmonic force constant FB is fitted',
     * ' using the second',/,1X,T5,'derivative w.r.t. the angle at the',
     * ' equilibrium collinear geometry.',/,1X,T5,'If FB is',
     * ' positive the quartic force constant is fitted using the',
     * /,1X,T5,'fourth derivative w.r.t. the angle.  If FB is',
     * ' negative, the harmonic',/,1X,T5,'and quartic force constants ',
     * 'are fitted to reproduce the location and ',
     * /,1X,T5,'depth of the well.')
624   FORMAT (1X,T5,'IBOPT=3, the harmonic and quartic bend force',
     * ' constants are fit to the',/,1X,T5,'potential at two angles;',
     * ' THETA=', F10.5,' and ',F10.5)
626   FORMAT (1X,T5,'IGS(6) = 2, the semiclassical approximation ',
     * 'is used',/,1X,T5,'to determine the eigenvalues of the ',
     * 'harmonic-quartic bending potential.')
628   FORMAT (1X,T5,'IGS(6) = 3, for FB<0, the eigenvalues of the',
     * ' harmonic-quartic potential',/,1X,T5,'are determined by the ',
     * 'WKB approximation.',/,1X,T5,'For FB>0, the eigenvalues of ',
     * 'the harmonic-quartic potential are ',
     * /,1X,T5,'determined by the perturbation-variation method.')
630   FORMAT (1X,T5,'Harmonic approximation used for the bend.')
632   FORMAT (/,1X,T2,'The state-selected option has been chosen ',
     * 'for the stretch.',/,1X,T2,'Rate constants are state-selected ',
     * 'ones for NSTR = ', I5)
634   FORMAT (/,1X,T2,'The state-selected option has been chosen for ',
     * 'the bends.',/,1X,T2,'Rate constants are state-selected ones ',
     * 'for NBEND = ', I5, ' and ', I5)
636   FORMAT (/,1X,T2,'Options for the LAG calucation:')
638   FORMAT (1X,T5,'The LAG calculations are state specific, summing',
     * ' over',/,1X,T5,'final vibrational states from', I4, ' to', I4)
640   FORMAT (1X,T5,'The period of the stretching vibration is',
     * ' approximated using the',1X,T5,'zero-point energy level.')
642   FORMAT (1X,T5,'The period of the stretching vibration is',
     * ' calculated using the',
     * /,1X,T5,'derivative of the phase integral.')
646   FORMAT (1X,T5,'The MEP path is the zero-order one.')
648   FORMAT (1X,T5,'For the LAG integrals over s, segments of size',
     *        F10.6,/,1X,T5,'are used with ', I3, 
     *              '-order Gauss quadrature in each segment.',
     *      /,1X,T5,'The order of the interpolation is ',I3)
650   FORMAT (1X,T5,'A root search is not performed for the ',
     *              'optimum alpha value, but',
     *      /,1X,T5,'it is chosen from the values 0.0, ',
     *              '0.5, and 1.0.')
652   FORMAT (1X,T5,'In the search for the alpha values that ',
     *              'minimize the phase integral',
     *      /,1X,T5,'the error tolerance in P(E) is ', 1PE13.5)
654   FORMAT (/,1X,35('*'),1X,'Masses',1X,35('*'),
     *       //,1X,T5, 'Atomic masses (atomic units)', T40, 
     *                 'Atom', T52,  A2, T60, F18.9,
     *      /,1X,T52, A2, T60, F18.9, /, 1X, T52, A2, T60,F18.9)
656   FORMAT (/,1X, T5, 'Reduced masses (atomic units)', T40, 2A2, 
     *                  ' pair (M) ', T60, F18.9,
     *        /,1X,T40,2A2,' pair (MP)', T60, F18.9,
     *        /,1X,T40,A2,',',2A2,' pair (MU)', T60, F18.9,
     *      //,1X, T5, 'Scaling factors (unitless)', T40, 'C/(B+C)', 
     *             T60, F18.9, 
     *       /,1X, T40, 'sqrt(M/MU)', T60, F18.9,
     *       /,1X, T5, 'Skew angle (degrees)', T60, F18.9)
658   FORMAT (/,1X,78('*'),//,1X,T5, 'Energy at bottom of asymptotic ',
     *                               'BC well defines', T55, 'zero',
     *        /,1X,T5, 'Energy of three atoms infinitely separated', 
     *             T50, 1PE17.8, T70, 'Hartree', /, T50, 1PE17.8, 
     *             T70, 'eV', /, T50, 1PE17.8, T70, 'kcal',
     *        /,1X,T5, 'Endoergicity', T50, 1PE17.8, T70, 'Hartree', 
     *        /,1X,T50, 1PE17.8, T70, 'eV', 
     *        /,1X,T50, 1PE17.8, T70, 'kcal')
660   FORMAT (/,1X,T5, 'Warning: There is no saddle point for ',
     *                 'this potential')
661   FORMAT (/,1X,T5,'IGS(2)=4 option selected.  The contribution of ',
     *        'the electronic partition function',/,1X,T5,'will not ',
     *        'be included in the 1D rate constant.')
662   FORMAT (/,1X,28('*'),2X, 'Partition functions', 1X, 28('*'),
     *       //,1X,T5, 'QTR - translational partition function in ',
     *                 'units of molec/cm (collinear)',
     *        /,1X,T11, 'or molec/cc (3D)', /, 1X, T16, 'R',
     *        /,1X,T5, 'Capital Phi (T) - total 3D reactant ',
     *                 'partition function in units of molec/cc',
     *       //,1X,T5, A9, ' Vibrational partition functions used ',
     *                     'for stretch',
     *        /,1X,T5, 'Rotational partition functions include ',
     *                 'the rotational symmetry factors')
664   FORMAT (1X,T5, 'Stretch vibrational partition functions are ',
     *               'for the state ',
     *      /,1X,T5, 'selected reaction with NSTR = ', I5)
666   FORMAT (/,T28, 'ln', T40, 'ln', T51, 'ln', T77, 'ln', T97,
     * 'ln', T101, 'R', T117, 'ln', / 4X, 'Temp', 5X, 'Beta', T28,
     * 'QTR', T40, 'QTR', T51, 'QROT', T62, 'Rot', T67, 'Quantum',
     * T77, 'QVIB', T84, A9, T94, 'Cap.Phi (T)', T107,
     * 'Equilibrium constants',/,10X,'Hartree**-1',T25,'Collinear',
     * T40, '3D', T62, 'Sig', T68, 'Corr.', T87, 'Corr.', T97,
     * '(Morse)', T107, 'Collinear', T121, '3D'/)
668   FORMAT (1X, F8.2, 1PE12.4, 2X, 3E12.4, 0P,I4, F9.4, 1PE12.4, 
     *            0PF8.4, 1P,3E12.4)         
670   FORMAT (/,1X, 25('*'), 1X, 'Canonical conventional TST', 
     *          1X, 25('*'),
     *       //,1X, T5, A9,' approximation for stretching ',
     *                     'degree of freedom',
     *        /,1X, T5, A10, ' approximation for bending ',
     *                       'degree of freedom')
672   FORMAT (1X, T5, 'Rate constants are state selected ones ',
     *                'for NSTR = ', I5)
673   FORMAT (/,1X,T5, '1D rate constants (cm/molec-sec)')
674   FORMAT (/,1X,T5, '3D rate constants (cc/molec-sec)')
676   FORMAT (/,T33,' ln(partition functions)',/,T24,A9, T43,
     * A10, T65, 'Rot', T70, 'Quantal', T79, 'Wigner', T94, 'TST',
     * T104, 'TST/Wigner', /, T5, 'T(K)', T14, 'Stretch', T26, 'Corr.',
     * T34, 'Bend**2', T46, 'Corr.', T55, 'Rot', T65, 'Sig', T71,
     * 'Corr.', T80, 'Kappa', T93, 'K (T)', T106, 'K (T)'/)
678   FORMAT (1X, F8.2, 1X, 2(1PE13.5, 0PF7.3), 1PE13.5, 0P,I3, F9.3,
     *            G13.5, 1P,2E13.5)   
682   FORMAT (/,1X,T5, 'Error: IOPT(2), NSTATE = ', I5, 1X, I5, ',',
     *        /,1X,T12, 'but IOPT(2) must be > or = NSTATE')
684   FORMAT (/,1X,T5, 'Error: IOPT(2) and IOPT(3) = ', I5, 
     *                 ', ', I5, ',',
     *        /,1X,T5, 'but IOPT(3) must be > or = IOPT(2)')
6001  FORMAT (/,1X,T5,'Error: Illegal input value for IGS(',I2,')')
6002  FORMAT (/,1X,T5,'Error: The value of the potential in the ',
     *                'asymptotic reactant region ',
     *       /,1X,T12,'is not zero, but is',1PE20.10,
     *       /,1X,T12,'The potential routine should be altered to ',
     *                'give zero in this region')
6003  FORMAT (/,1X,T5,'Error: For a system with no saddle point the ',
     *                'reaction must be run in the ',
     *        /,1X,T12,'exoergic direction')
6004  FORMAT (/,1X,T5,'Error: There is a problem with VMIN in the ',
     *                'call to set up the asymptotic ',
     *        /,1X,T12,'value to use to start reaction coordinate.')
      END
