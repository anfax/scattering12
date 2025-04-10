!!!*************************************************************
! 文件/File: wkbset.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: wkbset.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE WKBSET (IS, S, NLVL, E, UL, UG, PER)
C
C     WKBSET - set up grid of WKB energy levels
C        and entry points
C     WKBINT - interpolate WKB energy levels from grid
C     WKBRD  - read in grid of WKB energy levels from file 14
C     WKBWRT - write out grid of WKB energy levels to file 14
C
C  Called by:
C     ADIAB  - compute adiabatic potential curves
C     DERS   - derivatives of morse turning point and zero pt. energy
C              w.r.t. s
C     ESTR   - compute stretch vibrational energy levels
C     GEOM   - find equilibrium geometry for reactants and products
C     PFCNST - compute stretching part. fcn.
C     QMPART - compute part. fcn. for stretch
C     RPHSUM - summarize reaction path info
C     SADDLE - find saddle point(s) and do normal mode analysis
C     SAGCMP - compute info needed for effective mass terms
C     VIBTAU - compute vibrational period for u motion
C
C  Calls:
C     AITKF2 - Aitken interpolation
C     LOCS   - locate position of s in grid
C     WKB    - compute WKB energy levels for stretch
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LSET
      DIMENSION F(7)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      PARAMETER (NWKBDM=NSDM+6)
      DIMENSION ELVL0(NWKBDM), ELVLN(NWKBDM), ULTP(NWKBDM),
     *UGTP(NWKBDM), TAUS(NWKBDM)
      COMMON /INTER1/ NINT
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      CALL WKB (IS, S, NLVL, E, UL, UG,PER)
      INDEX = IS
      IF (IS .LT. 0) INDEX = NWKBDM + IS + 1
      IF (NLVL .NE. NSTATE) THEN
         ELVL0(INDEX) = E
         IF (LSYM .AND. IS .GE. 0) ELVL0(2*NMID-IS) = E
      ELSE
         ELVLN(INDEX) = E
         ULTP(INDEX) = UL
         UGTP(INDEX) = UG
         TAUS(INDEX) = PER
         IF (LSYM .AND. IS .GE. 0) THEN
            INDEX = 2*NMID - IS
            ELVLN(INDEX) = E
            ULTP(INDEX) = UL
            UGTP(INDEX) = UG
            TAUS(INDEX) = PER
         END IF
      END IF
      RETURN
C
      ENTRY WKBINT (IS, S, NLVL, E, UL, UG, PER)
C
      LSET = .FALSE.
      IF (IS .EQ. 0 .OR. IS .GT. NSMAX) THEN
         IF (S .LE. SS(1)) THEN
            INDEX = 1
         ELSE IF (S .GT. SS(NSMAX)) THEN
            INDEX = NSMAX
         ELSE
            LSET = .TRUE.
            CALL LOCS (INDEX, S)
            NP = NINT+1
            NN = NP/2 - 1
            INDEX = MAX(1, INDEX-NN)
            IF (INDEX+NP .GT. NSMAX) INDEX = NSMAX - NP
            IF (NLVL .NE. NSTATE) THEN
               E = AITKF2(S, ELVL0(INDEX), F, SS(INDEX), NINT)
            ELSE
               E = AITKF2(S, ELVLN(INDEX), F, SS(INDEX), NINT)
               UL = AITKF2(S, ULTP(INDEX), F, SS(INDEX), NINT)
               UG = AITKF2(S, UGTP(INDEX), F, SS(INDEX), NINT)
               PER = AITKF2(S, TAUS(INDEX), F, SS(INDEX), NINT)
            END IF
         END IF
      ELSE
         INDEX = IS
         IF (IS .LT. 0) INDEX = NWKBDM + 1 + IS
      END IF
      IF (.NOT.LSET) THEN
         IF (NLVL .NE. NSTATE) THEN
            E = ELVL0(INDEX)
         ELSE
            E = ELVLN(INDEX)
            UL = ULTP(INDEX)
            UG = UGTP(INDEX)
            PER = TAUS(INDEX)
         END IF
      END IF
      RETURN
C
      ENTRY WKBWRT
C  Write out WKB energy levels, etc. to unit 14
      OPEN (UNIT=14, FILE='abc.14', FORM='FORMATTED',                   GCL1096
     *      STATUS='UNKNOWN', ERR=910)                                  GCL1096
      WRITE (14, 1400) NWKBDM
      WRITE (14, 1401) (ELVL0(I), ELVLN(I), ULTP(I), UGTP(I), TAUS(I),
     *   I=1,NWKBDM)
      RETURN
C
910   WRITE (6, 911)                                                    GCL1096
      STOP 'WKBWRT 1'
911   FORMAT(/,1X,T5,'Error: Cannot open abc.14 in WKBWRT')             GCL1096
C
      ENTRY WKBRD
C  Read in WKB energy levels, etc. from unit 14
      OPEN (UNIT=14, FILE='abc.14', FORM='FORMATTED',                   GCL1096
     *      STATUS='OLD', ERR=900)                                      GCL1096
      READ (14, 1400) NNN
      IF (NNN.GT.NWKBDM) THEN
         WRITE (6, 6000)
         STOP 'WKBRD 1'
      END IF
      READ (14, 1401) (ELVL0(I), ELVLN(I), ULTP(I), UGTP(I), TAUS(I),
     *   I=1,NNN)
      RETURN
C
900   WRITE (6, 901)                                                    GCL1096
      STOP 'WKBRD 2'                                                    GCL1096
901   FORMAT(/,1X,T5,'Error: Cannot open abc.14 in WKBRD')              GCL1096
 1400 FORMAT (1X, I5)
 1401 FORMAT (1X, G19.10, 3G20.10)
6000  FORMAT(/,1X,T5,'Error: Number of energy levels exceeds the ',
     *               'dimensions of the arrays')
      END
