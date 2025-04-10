!!!*************************************************************
! 文件/File: vtstm.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vtstm.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      PROGRAM VTSTM
C
C     VTSTM  - main program
C
C  Calls:
C     ADIAB  - compute adiabatic potential curves
C     DATAIN - read in data and set up constants
C     GTST   - compute free energies, CVT and ICVT rates
C     KAPVA  - compute kappas
C     RPATH  - follow reaction path and store MEP info on grid
C     PORCPU - returns cpu time
C     SUMMRY - summarize rate constants
C     TIMING - prints out timing information, current and elasped time.
C
C  A potential routine must be supplied of the form as follows:
C     SUBROUTINE PREPEF (no arguments): must be supplied which is called
C        once to set up necessary parameters
C     SUBROUTINE PEF (R1, R2, R3, V, D1, D2, D3, NDER)
C           R1 - AB distance
C           R2 - BC distance
C           R3 = R1 + R2
C           A + BC - reactants
C           AB + C - products
C           V  - Energy
C           D1 - DV/DR1
C           D2 - DV/DR2
C           D3 - DV/DR3
C           NDER - integer flag that controls the derivative calculations
C                = 0, no derivatives are needed
C                = 1, first derivatives are needed
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NQGKDM=81)
      PARAMETER (NSADDM=4)
      COMMON /ESPEC/  ESPEC(40), NESPEC
      LOGICAL LGS(10), LGS2(10), LLAG, LLAGRS, LRSTRT
      COMMON /LOGIC/  LGS
      COMMON /LOGIC2/ LGS2
      COMMON /LAGCOM/ PT3(NQGKDM), WT3(NQGKDM,2), NQ32, NSEG3, IOPTAU,
     * LLAG, LLAGRS
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     * NSSP(NSADDM), LRSTRT
      CHARACTER*20 CFSTAT
      SAVE CFSTAT                                                       TCA1097
      DATA CFSTAT /'KEEP'/
C
C   Open the input and output files
C
      OPEN (UNIT=5,FILE='abc.5',FORM='FORMATTED',STATUS='OLD',ERR=10)   GCL1096
      OPEN (UNIT=6,FILE='abc.6',FORM='FORMATTED',                       GCL1096
     *      STATUS='UNKNOWN',ERR=20)                                    GCL1096
C
      CALL HEADR                                                        GCL1092
      CALL PORCPU(TIMEI)                                                GCL0992
      WRITE (6, 600) TIMEI
C
C  DATAIN reads in data, sets up constants, and does conventional TST
C     calculation.  All inputs are through DATAIN except parameters for
C     the potential subroutine which may be read in through PREPEF.  All
C     units are internally in atomic units except where otherwise noted.
C
C  Option variables in COMMON /LOGIC/ ,/LOGIC2/, and /OPTION/ all are
C     read in thru DATAIN:
C     LGS(1) - if .true. print is on in datain
C        (2) - if .true. execution is terminated after datain
C        (3) - if .true. morse correction to stretching vibration
C              is used.
C        (4) - if .true. quartic correction to bending vibration
C              used
C        (5) - if .true. execution is terminated after the
C              vibrationally adiabatic potential curve is computed
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
C    IOPT(6-20) - not used
C
      CALL DATAIN
C
C Close files 
      CFSTAT = 'DELETE'                                                 GCL1096
      CLOSE (UNIT=5, STATUS=CFSTAT)                                     GCL1096
      IF (LGS2(5)) CLOSE (UNIT=14, STATUS=CFSTAT)                       GCL1096
      CFSTAT = 'KEEP'                                                   GCL1096
      CALL TIMING (TIMEF, TIMEI)
C
      IEND = 0
      IF (LGS(2)) THEN
         IEND = 1
      ELSE
C  Find reaction coordinate by starting at the saddle point and
C     following the gradient down.  The initial step is taken along the
C     normal coordinate at the saddle point which has an imaginary
C     frequency.
         CALL RPATH
         CALL TIMING (TIMEF, TIMEI)
C  Compute adiabatic potential curves
            CALL ADIAB
C Close files 
            IF (LGS2(7)) CLOSE (UNIT=1, STATUS=CFSTAT)                  GCL1096
            IF (LRSTRT) THEN                                            GCL1096
                CFSTAT = 'DELETE'                                       GCL1096
                CLOSE (UNIT=1, STATUS=CFSTAT)                           GCL1096
                CFSTAT = 'KEEP'                                         GCL1096
            ENDIF                                                       GCL1096
            IF (LGS2(4)) CLOSE (UNIT=14, STATUS=CFSTAT)                 GCL1096
            CALL TIMING (TIMEF, TIMEI)
C
         IF (LGS(5)) THEN
            IEND = 2
         ELSE
C  Compute transmission coefficients
            CALL KAPVA
C Close files 
            IF (LLAG .AND. LLAGRS) THEN                                 GCL1096
                CFSTAT = 'DELETE'                                       GCL1096
                IF (LGS(9)) CLOSE (UNIT=11, STATUS=CFSTAT)              GCL1096
                IF (LGS(10)) CLOSE (UNIT=12, STATUS=CFSTAT)             GCL1096
                CFSTAT = 'KEEP'                                         GCL1096
            ENDIF                                                       GCL1096
            IF (LLAG .AND. LGS2(8)) THEN                                GCL1096
                IF (LGS(9)) CLOSE (UNIT=11, STATUS=CFSTAT)              GCL1096
                IF (LGS(10)) CLOSE (UNIT=12, STATUS=CFSTAT)             GCL1096
            ENDIF                                                       GCL1096
            CALL TIMING (TIMEF, TIMEI)
C
            IF (LGS(6)) THEN
               IEND = 3
            ELSE
               IF (NESPEC .LE. 0) THEN
C  Free energy calculation (canonical variational tst)
                  CALL GTST
                  CALL TIMING (TIMEF, TIMEI)
C  Summarize results
                  CALL SUMMRY
                  CALL TIMING (TIMEF, TIMEI)
               END IF
            END IF
         END IF
      END IF
      IF (IEND .NE. 0) WRITE (6,602) IEND
C
      STOP
C
10    CONTINUE                                                          GCL1096
      WRITE (6, 900)                                                    GCL1096
      STOP 'VTSTM 1'                                                    GCL1096
20    CONTINUE                                                          GCL1096
      WRITE (6, 901)                                                    GCL1096
      STOP 'VTSTM 2'                                                    GCL1096
C
  600 FORMAT (/,1X, 78('*'),/,1X,T5,'Start the calculation, ',
     *   'the time is', F12.2, ' seconds ',/,1X,78('*'))
  602 FORMAT (/,1X, 78('*'),/,1X,T5,'IGS(1) was set equal to ',I1,      TCA1296
     *        ', the calculation is now terminated',/,1X,78('*'))
900   FORMAT(/,1X,'Error: Cannot open input file abc.5')                GCL1096 
901   FORMAT(/,1X,'Error: Cannot open output file abc.6')               GCL1096
C
      END
C
