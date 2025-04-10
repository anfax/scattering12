!!!*************************************************************
! 文件/File: tabl21.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: tabl21.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE TABL21 (T, IS, ITYP, IC3D)
C
C     TABL21 - print out table of GTS info
C
C     Modified 3/18/91 to include centrifugal oscillator bend energies
C
C  Called by:
C     ADIAB  - compute adiabatic potential curves
C     EXTRAS - print out extra info about GTS
C
C  Calls:
C     COBEND - compute semiclassical eigenvalue of centrifugal
C              oscillator
C     EBEND  - compute bending energy levels
C     TITLES - print 2-line title
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LFIRST, LWRITE, LERR
      CHARACTER*3 ATYP(2)
      CHARACTER*2 AA(2)
      CHARACTER*1 AVA
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      LOGICAL LGS2(10)
      COMMON /LOGIC2/ LGS2
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /MORAB/  DAB, XKAB, AMAB, VDELTA
      COMMON /MORBC/  DBC, XKBC, AMBC
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      COMMON /STATE2/ DGBND, LSBEND, NBND1, NBND2
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
      SAVE LFIRST, LWRITE, LERR, AA, ATYP                               TCA1097
      DATA LFIRST, LWRITE /.TRUE., .TRUE./
      DATA LERR /.FALSE./                                               GCL94
      DATA AA /'1D', '3D'/
      DATA ATYP/'a  ','CVT'/
C
      IF (.NOT.LWRITE .OR. .NOT.LGS2(1)) RETURN
      IF (LFIRST) THEN
          OPEN (UNIT = 21, FILE='abc.21', FORM='FORMATTED',             GCL1096
     *          STATUS='NEW', ERR=100)                                  GCL1096
         GO TO 20
   10    CONTINUE
            LERR = .TRUE.
   20    CONTINUE
         IF (LERR) THEN
            WRITE (6, 6000)
            LWRITE = .FALSE.
            RETURN
         END IF
         CALL TITLES (1, 21, 1)
         IF (LSTATE .NE. 0) WRITE (21, 2105) NSTATE
         AVA = ' '
         IF (LSTATE .EQ. 0 .OR. NSTATE.EQ.0) AVA = 'G'
         IF (LSTATE .NE. 0 .AND. NSTATE .GT. 0 .AND. NSTATE .LT. 10)
     *      AVA = CHAR(48+NSTATE)
         ATYP(1)(2:2) = AVA
         WRITE (21, 2101) AVA, AVA
         LFIRST = .FALSE.
      END IF
      IF (IS .LE. 0) THEN
         S = -50.0D0
         R1 = 50.0D0
         R2 = R2ASY
         VMEP = 0.D0                                                    GCL1092
         D = DBC
         XK = XKBC
         WB = 0.D0                                                      GCL1092
         EB = 0.D0                                                      GCL1092
      ELSE IF (IS .GT. NSDM10) THEN
         S = 50.D0                                                      GCL1092
         R1 = R1ASY
         R2 = 50.D0                                                     GCL1092
         VMEP = VDELTA*CKCAL
         D = DAB
         XK = XKAB
         WB = 0.D0                                                      GCL1092
         EB = 0.D0                                                      GCL1092
      ELSE
         S = SS(IS)
         D = DS(IS)
         XK = XKS(IS)
         R2 = Y2(IS)/CM2
         R1 = X2(IS) - CM1*R2
         VMEP = VS(IS)*CKCAL
         WB = 0.D0                                                      GCL1092
         IF (FBS(IS) .GT. 0.D0) WB = SQRT(FBS(IS)*GBS(IS))*CCM
         IF (LGS2(6)) THEN
C  centrifugal oscillator energy level for bend
            CALL COBEND(NBND1,NBND2,FBS(IS),QFBS(IS),GBS(IS),EB)
         ELSE
C  uncoupled bending energy level
            EB = EBEND(NBND1, FBS(IS), QFBS(IS), GBS(IS))
            IF (NBND1 .EQ. NBND2) EB = 2.D0*EB
            IF (NBND1 .NE. NBND2) EB = EB + EBN(NBND2, FBS(IS),
     *          QFBS(IS),GBS(IS))
         END IF
         EB = EB*CKCAL
      END IF
      WSTR = 4.D0*D*CCM/XK
      ES = ESTR(IS, S, D, XK)*CKCAL
      VA1D = VMEP + ES
      VA3D = VA1D + EB
      IF (ITYP .EQ. 0) THEN
         WRITE (21, 2104) S, R1, R2, VMEP, VA1D, VA3D, WSTR, ES, WB, EB
      ELSE IF (ITYP . EQ. 1) THEN
         WRITE (21, 2102) ATYP(1), AA(IC3D), S, R1, R2, VMEP, VA1D,
     *      VA3D, WSTR, ES, WB, EB
      ELSE IF (ITYP .EQ. 2) THEN
         WRITE (21, 2103) ATYP(2), AA(IC3D), S, T, R1, R2, VMEP, VA1D,
     *      VA3D, WSTR, ES, WB, EB
      END IF
      RETURN
C
100   WRITE (6, 900)                                                    GCL1096
      LWRITE = .FALSE.                                                  GCL1096
C
 2101 FORMAT (/,1X,T5,'Bottleneck properties',
     *        /,1X,T5,'Distances in bohr, energies in kcal, and ',
     *                'frequencies in cm**-1',//, T26, 'T(K)',
     *   T37, 'R1', T44, 'R2', T53, 'V', 4X, 'Va', A1, '(1D) Va', A1,
     *   '(3D)', T77, 'WSTR', T84, 'ESTR', T92, 'WBND', T99, 'EBND'/)
 2102 FORMAT (2X, 's*', A3, ' (', A2, ')=', F7.3, 12X, 2F7.2, 2X,
     *   3(F7.2, 1X), 2(2X, F6.0, F7.2))
 2103 FORMAT (2X, 's*', A3, ' (', A2, ')=', F7.3, 2X, F8.2, 2X, 2F7.2,
     *   2X, 3(F7.2, 1X), 2(2X, F6.0, F7.2))
 2104 FORMAT (2X,'s*TST', 5X, '=', F7.3, 12X, 2F7.2, 2X, 3(F7.2, 1X), 
     *        2(2X,F6.0, F7.2))
 2105 FORMAT (/, 1X, T5, 'State-selected option is chosen for ',
     *                   'NSTR = ', I3, '.')
900   FORMAT(/,1X,T5,'Error: Cannot open abc.21 in TABL21')
6000  FORMAT(/,1X,T5,'Error: In TABL21 there is a problem writing ',
     *               'the bottleneck properties ',
     *       /,1X,T5,'to the file linked to FORTRAN unit 21')
      END
