!!!*************************************************************
! 文件/File: gtst.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: gtst.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE GTST
C
C     GTST   - compute free energies, CVT, ICVT AND CUS rates
C
C     Modified 3/18/91 to include centrifugal oscillator bend energies
C
C  Called by:
C     VTSTM  - main program
C
C  Calls:
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C     EBEND  - compute bending energy levels
C     EXTRAS - print out extra info about GTS
C     EXTREM - finds extrema on grid and stores locations in grid
C     INTERP - interpolate r.p. info from grid
C     PFCN   - compute GTS partition functions and free energy
C     QUADFT - quadratic fit of three points
C     STORE  - store reaction path info
C     THRCOR - compute threshold corrections for ICVT
C     TITLES - print out 2-line title
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LIND(2), LCONV, LPRNT
      CHARACTER*9 AA(3)
      CHARACTER*15 AH(3), AB
      DIMENSION B(3), GEXTR(4, 2), GX(2), GS(3, 4), GMAX(4), GSP(2),
     *SGN(2), QTOT(2), VMXCVT(2)
      DIMENSION ISEXTR(4, 4), NMAX(2), NEXTR(2), ISMX(4)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /ADIAB1/ VAD(NSDM,2), VMAX(2), VR, VP, SMAX(2), ISMAX(2)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /C12/    C12, CONC12, INDFR
      COMMON /EACT1/   IACT, NT1(10), NT2(10)
      PARAMETER (NTEMDM=100)                                            TCA0996
      DOUBLE PRECISION KAPW
      LOGICAL PTEMP                                                     TCA0996
      COMMON /TEMPCM/ TEMP(NTEMDM), BETA(NTEMDM), CPHI(NTEMDM),
     *CNST(NTEMDM), CPHIC(NTEMDM), CNSTC(NTEMDM), EQUIL(NTEMDM,2),
     *KAPW(NTEMDM,2), RATIO(NTEMDM,2), NTMAX, PTEMP(NTEMDM)             TCA0996
      COMMON /ELFCT/  ELFACT(NTEMDM,2)
      COMMON /GTST1/   NPGTST
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      LOGICAL LGS2(10)
      COMMON /LOGIC2/  LGS2
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /MORAB/  DAB, XKAB, AMAB, VDELTA
      COMMON /MORBC/  DBC, XKBC, AMBC
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /MORSTS/ DTS(NSADDM), XKTS(NSADDM), AMTS(NSADDM),
     *OMTS(NSADDM), OMIMG(NSADDM)
      COMMON /OPTION/ IOPT(20)
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
      COMMON /PFCN1/  QROT, CROT, QB, CB, DELPHI, CSTR, QSTR, QT, G(2)
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     *XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     *NSAD, NSADMX(2)
      COMMON /STATE2/ DGBND, LSBEND, NBND1, NBND2
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      PARAMETER (NKATYP=21)
      PARAMETER (NOB=4)
      PARAMETER (NKALAB=NKATYP-NOB)
      PARAMETER (NKTYP=NKALAB+1)
      PARAMETER (NRATE=5)
      CHARACTER*10 ATYP(NKALAB)
      LOGICAL LPTYP(NKALAB)                                             TCA0197
      COMMON /TYPES/  LPTYP, ATYP                                       TCA0197
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C  KAP(IT,J,I), IT - index over temperatures, J - index over methods,
C     I - 1=1D, 2=3D
      DOUBLE PRECISION KAP
      COMMON /BLKKCM/ KAP(NTEMDM,NKATYP,2)                              GCL1096
      COMMON /GTSTCM/ XK(NTEMDM,NRATE,2), XKPOL(NTEMDM,2),              GCL1096
     *                QSV(NTEMDM,4,5), THCOR(NTEMDM,2), RAT(4)          GCL1096 
      DIMENSION GSV(NSDM10,2), QTSV(NSDM10,2)                           GCL1096
      EQUIVALENCE (QTOT, QSTR)
      SAVE AVAG, C6, AA, AB, AH                                         TCA1097
      DATA AVAG /6.022169D23/, C6 /0.16666666666667D0/
      DATA AA /'Harmonic ', '  Morse  ', 'WKB-Morse'/, AB /' '/
      DATA AH /'Saddle point', 'Collinear', '3D'/
C
      IF (LGS(3)) THEN
         IF (LGS(7)) THEN
            IMH = 3
         ELSE
            IMH = 2
         END IF
      ELSE
         IMH = 1
      END IF
      WRITE (6, 600)
      CALL TITLES (1, 6, 1)
      WRITE (6, 602) AA(IMH)
      IF (LGS(4)) THEN
         WRITE (6, 604)
      ELSE
         WRITE (6, 605)
      END IF
      N = NSADMX(2)
      IF (LGS(9)) THEN
         IC3D1 = 1
      ELSE
         IC3D1 = 2
      END IF
      IF (LGS(10)) THEN
         IC3D2 = 2
      ELSE
         IC3D2 = 1
      END IF
      IF (LSYM) THEN
         NSLAST = NMID
      ELSE
         NSLAST = NSMAX
      END IF
      IF (NSAD .GT. 0) THEN
         N1 = NSSP(NSADMX(1))
         N2 = NSSP(NSADMX(2))
      END IF
C  Compute rate constants with dividing surface at maximum of
C    ground-state adiabatic potential (AG rates)
      DO 400 IT = 1,NTMAX
         BET = BETA(IT)
         RT = CKCAL/BET
         T = -LOG(BET*TAU*TPI)
         DO 400 IC3D = IC3D1,IC3D2
            CALL PFCN(ISMAX(IC3D),IT,AB,.FALSE.)
            XK(IT,5,IC3D) = EXP(-G(IC3D)/RT+T)
  400 CONTINUE
C  Loop over temperatures
      DO 300 IT = 1,NTMAX
         BET = BETA(IT)
         RT = CKCAL/BET
         IF (PTEMP(IT)) THEN                                            TCA0996
           WRITE (6, 606) TEMP(IT), BET, AA(IMH)
         ELSE                                                           TCA0996
           WRITE (6, 630) TEMP(IT)                                      TCA0996
         END IF                                                         TCA0996
C  free energy at reactants
         IF (IOPT(4) .NE. 0 .AND. PTEMP(IT)) WRITE (6, 607)             TCA0996
         CALL PFCN (0, IT, AB, .TRUE.)
         IF (IOPT(4) .EQ. 0) THEN
            GX(1) = G(1)
         ELSE
            GX(1) = -1.D35
         END IF
         GX(2) = -1.D35
C   search for extrema in G(T,S).  GEXTR(I,J) contains extrema J for
C   col(i=1) and 3d(i=2).  J=1,3 are max., J=2 is min
         LIND(1) = .FALSE.
         LIND(2) = .FALSE.
         IPRNT = NPGTST
         ISAD = 0
         ISP = 0
C  loop over grid points
         DO 20 IS = 1,NSLAST
C  check if information is to be printed
            IF (IS .GE. ISP) THEN
               ISAD = ISAD + 1
               IF (ISAD .LE. NSAD) THEN
                  ISP = MIN (NSLAST, NSSP(ISAD))
               ELSE
                  ISP = NSLAST
               END IF
               LPRNT = .TRUE.
            ELSE
               IPRNT = IPRNT + 1
               LPRNT = IPRNT .GE. NPGTST
            END IF
            IF (LPRNT) IPRNT = 0
            CALL PFCN (IS, IT, AB, LPRNT)
C  check if free energy is an extremum
            DO 10 I = IC3D1,IC3D2
               GSV(IS,I) = G(I)
               QTSV(IS,I) = QTOT(I)
               CALL EXTREM (IS, GX(I), G(I), LIND(I), SGN(I), NMAX(I),
     *            NEXTR(I), ISEXTR(1,I), GEXTR(1,I))
               GX(I) = G(I)
   10       CONTINUE
            IF (NSAD .GT. 0) THEN
               IF (IS .EQ. N1) GSP(1) = G(1)
               IF (IS .EQ. N2) GSP(2) = G(2)
            END IF
   20    CONTINUE
C  free energy in product region
         IF (LSYM) THEN
C  symmetric system
            IF (PTEMP(IT)) WRITE (6, 608)                               TCA0996
            IS = NSLAST + 1
            DO 30 I = IC3D1,IC3D2
               G(I) = GSV(NMID-1,I)
               CALL EXTREM (IS, GX(I), G(I), LIND(I), SGN(I), NMAX(I),
     *            NEXTR(I), ISEXTR(1,I), GEXTR(1,I))
               IF (NEXTR(I) .EQ. 1) GO TO 30
               IF (NEXTR(I) .GT. 2 .AND. GEXTR(1,I) .LT. GEXTR(3,I))
     *            THEN
                  N = 3
               ELSE
                  N = 1
               END IF
C  check if abs max is at sym point
               IF (ISEXTR(N,I) .EQ. NMID) GO TO 30
               IF (N .EQ. 1) THEN
C  abs max is at first location (N=1)
                  GEXTR(3,I) = GEXTR(1,I)
                  ISEXTR(3,I) = 2*NMID - ISEXTR(1,I)
                  IF (NEXTR(I) .LE. 3 .OR. GEXTR(2,I) .LT. GEXTR(4,I))
     *               GO TO 30
                  GEXTR(2,I) = GEXTR(4,I)
                  ISEXTR(2,I) = ISEXTR(4,I)
               ELSE
C  abs max is at third location (N=3)
                  GEXTR(1,I) = GEXTR(3,I)
                  ISEXTR(1,I) = ISEXTR(3,I)
                  ISEXTR(3,I) = 2*NMID - ISEXTR(3,I)
                  GEXTR(2,I) = GEXTR(4,I)
                  ISEXTR(2,I) = ISEXTR(4,I)
               END IF
   30       CONTINUE
         ELSE
C  nonsymmetric system
            CALL PFCN (NSDM10+1, IT, AB, .TRUE.)
            DO 40 I = IC3D1,IC3D2
               IS = NSMAX + 1
               IF (IOPT(5) .EQ. 0) THEN
                  CALL EXTREM (IS, GX(I), G(I), LIND(I), SGN(I),
     *               NMAX(I), NEXTR(I), ISEXTR(1,I), GEXTR(1,I))
                  IS = NSMAX + 2
                  GX(I) = G(I)
               ELSE
                  IF (PTEMP(IT)) WRITE (6, 609)                         TCA0996
               END IF
               G(I) = -1.D35
               CALL EXTREM (IS, GX(I), G(I), LIND(I), SGN(I), NMAX(I),
     *            NEXTR(I), ISEXTR(1,I), GEXTR(1,I))
   40       CONTINUE
         END IF
C  threshold correction (ICVT)
         DO 90 I = IC3D1,IC3D2
            IP = I + 2
            NN = MIN(NEXTR(I), 3)
C  loop over maxima in the uncorrected free energy curves
            DO 80 IEX = 1,NN,2
               IS = ISEXTR(IEX, I)
               IF (IS .GT. NSMAX) THEN
                  IS = NSDM10+1
                  ISEXTR(IEX,I) = IS
               END IF
               ISM = IS
C  skip correction for a maximum off the grid, at the left most edge of
C     the grid, or at the right most edge of the grid of a nonsymmetric
C     system.
               IF (IS .LE. 1 .OR. IS .EQ. NSMAX .OR.
     *            IS .GT. NSDM10) GO TO 70
C  skip correction for a maximum past the middle of the grid of a
C     symmetric system.
               IF (LSYM .AND. IS. GT. NMID) GO TO 70
C  compute correction for three points centered around uncorrected
C     maximum
               IS = IS - 2
               GM = -1.D35
               LCONV = .TRUE.
               DO 50 J = 1,3
                  IS = IS + 1
C  if IS>NMID for a symmetric system it must be for J=3 and will be the
C     same as for J=1, therefore skip it
                  IF (LSYM .AND. IS .GT. NMID) GO TO 50
                  CALL THRCOR (I, IS, BET, TCOR, QTSV(IS,I))
                  IF (TCOR .LE. 0.0D0) GO TO 70
                  GT = GSV(IS,I) - RT*LOG(TCOR)
C  check for maximum on the grid of 3 points
                  IF (GT .GE. GM) THEN
                     ISM = IS
                     JMAX = J
                     GM = GT
                  END IF
                  IF (J .GT. 1) THEN
                     T = GT - G2
                     IF (ABS(GT) .GT. 1.D-8) T = T/GT
                     LCONV = LCONV .AND. ABS(T) .LT. 1.D-8
                  END IF
                  G2 = GT
   50          CONTINUE
C  If the maximum of the corrected free energy on the grid of 3 occurs
C  in the middle then we are close enough for the quadratic fit.  Also,
C  if the difference in free energy between successive points is small
C  then don't look any more.
               IF (JMAX .EQ.  2 .OR.  LCONV) GO TO 70
               IF (JMAX .EQ. 1) THEN
                  ISGN = -1
               ELSE
                  ISGN = 1
               END IF
               IS = ISM
C  search for maximum by skipping along the big grid in the correct
C  direction.
   60          CONTINUE
                  IF (IS .LE. 1 .OR. IS .GE. NSLAST) GO TO 70
                  IS = IS + ISGN
                  CALL THRCOR (I, IS, BET, TCOR, QTSV(IS,I))
                  IF (TCOR .LE. 0.0D0) GO TO 70
                  GT = GSV(IS,I) - RT*LOG(TCOR)
                  IF (GT .LT. GM) GO TO 70
                  T = GT - GM
                  IF (ABS(GT) .GT. 1.D-8) T = T/GT
                  IF (ABS(T) .LT. 1.D-8) GO TO 70
                  GM = GT
                  ISM = IS
               GO TO 60
   70          CONTINUE
               ISEXTR(IEX,IP) = ISM
   80       CONTINUE
   90    CONTINUE
C
C  Search for extrema.  Loop over 4 types (uncorrected 1D and 3D and
C     corrected 1D and 3D)
         DO 200 I = 1,4
            IF ((I .EQ. 1 .OR. I .EQ. 3) .AND. .NOT.LGS(9)) GO TO 200
            IF ((I .EQ. 2 .OR. I .EQ. 4) .AND. .NOT.LGS(10)) GO TO 200
            GMAX(I) = -1.D35
            IM2 = MOD(I-1, 2) + 1
            IF (I .LT. 3) THEN
               IF (PTEMP(IT)) WRITE (6, 610) AH(IM2+1)                  TCA0996
            ELSE
               IF (PTEMP(IT))WRITE (6, 612) AH(IM2+1)                   TCA0996
            END IF
            NN = MIN(3, NEXTR(IM2))
            NSTEP = 1
            IF (I .GT. 2) NSTEP = 2
C  Loop over extrema in curve.  For threshold corrected curves only
C     maxima.
            DO 150 IEX = 1,NN,NSTEP
               SGNG = (-1)**IEX
               ISM = ISEXTR(IEX,I)
               ISI = ISM - 1
C              WRITE (6, 6610) I, IEX, ISM, SS(ISM), SGNG
C6610          FORMAT (' GTST: I,IEX,ISM,S,SGNG=', 3I5, 1P2E15.7)
               IF (LSYM .AND. ISM .EQ. NMID) GO TO 140
               IF (LSYM .AND. ISM .GT. NMID) GO TO 160
C  Skip search if the extrema is at the edge or off the grid.
               IF (ISM .LE. 1 .OR. ISM .GE. NSLAST) GO TO 140
               GM = SGNG*1.D35
               ISM2 = ISM - 2
C  Compute free energy at 3 points centered at extrema
               LCONV = .TRUE.
               DO 100 J = 1,3
                  IS = ISM2 + J
                  GT = GSV(IS,IM2)
                  IF (I .GT. 2) THEN
                     CALL THRCOR (IM2, IS, BET, TCOR, QTSV(IS,IM2))
                     IF (TCOR .LE. 0.0D0) GO TO 140
                     GT = GT - RT*LOG(TCOR)
                  END IF
                  IF (SGNG*(GT-GM) .LE. 0.D0) THEN
                     JMAX = J
                     GM = GT
                     ISM = IS
                  END IF
C                 WRITE (6, 6611) J, JMAX, ISM, IS, SS(IS), GM, GT
C6611 FORMAT (' J,JMAX,ISM,IS,S,GM,GT=', 4I5, 1P3E15.7)
                  IF (J .GT. 1) THEN
                     T = GT - G2
                     IF (ABS(GT) .GT. 1.D-8) T = T/GT
                     LCONV = LCONV .AND. ABS(T) .LT. 1.D-8
                  END IF
                  G2 = GT
                  GS(J,I) = GT
  100          CONTINUE
C  search for maximum in curve, maximum of 10 iterations
               ISST = NSDM4
               ISSTD = 3
               IC = 0
               DELS = 0.5D0*(SS(ISI+1) - SS(ISI))
  110          CONTINUE
               IF (IC .GE. 10 .OR. LCONV) GO TO 130
                  IC = IC + 1
                  GQT = 0.D0                                            GCL1092
                  S = SS(ISM)
                  IF (JMAX .EQ. 2) THEN
C  extrema at mid point use quadratic fit
                     CALL QUADFT (SS(ISI), GS(1,I), B)
                     IF (SGNG*B(3) .GT. 0.0D0) THEN
                        S = -.5D0*B(2)/B(3)
                        IF (S .GT.  SS(ISI) .AND.  S .LT.  SS(ISI+2))
     *                     THEN
                           GQT = (B(2)+B(3)*S)*S + B(1)
                        ELSE
                           S = SS(ISM)
                        END IF
                     END IF
                  END IF
                  S = S - 2.0D0*DELS
                  JMAX = 0
C  compute free energy on new grid
                  LCONV = .TRUE.
                  DO 120 J = 1,3
                     S = S + DELS
                     CALL INTERP (ISI)
                     IS = ISST + J
                     CALL STORE (IS)
                     CALL PFCN (IS, IT, AB, .FALSE.)
                     GT = G(IM2)
                     IF (I .GT. 2) THEN
                        CALL THRCOR (IM2, IS, BET, TCOR, QTOT(IM2))
                        IF (TCOR .LE. 0.0D0) THEN
                           IS = ISM
                           GO TO 140
                        END IF
                        GT = GT - RT*LOG (TCOR)
                     END IF
                     IF (SGNG*(GT-GM) .LE. 0.D0) THEN
                        JMAX = J
                        ISM = IS
                        GM = GT
                     END IF
C                    WRITE (6, 6611) J, JMAX, ISM, IS, SS(IS), GM, GT
                     IF (J .GT. 1) THEN
                        T = GT - G2
                        IF (ABS(GT) .GT. 1.D-8) T = T/GT
                        LCONV = LCONV .AND. ABS(T) .LT. 1.D-6
                     END IF
                     G2 = GT
                     GS(J,I) = GT
  120             CONTINUE
                  IF (LCONV) GO TO 130
                  IF (JMAX .NE. 0) THEN
                     ISI = ISST + 1
                     ISST = ISST + ISSTD
                     ISSTD = -ISSTD
                     T = GM - GQT
                     IF (ABS(GM) .GT. 1.D-8) T = T/GM
                     LCONV = ABS(T) .LT. 1.D-6
                  END IF
                  IF (JMAX .EQ. 0 .OR. JMAX .EQ. 2) DELS =
     *               0.5D0*DELS
                  LCONV = LCONV .OR. DELS .LT. 1.D-9
               GO TO 110
  130          CONTINUE
C  extrema found
               S = SS(ISM)
  140          CONTINUE
C  compute free energy and print out info for extrema
               CALL PFCN (ISM, IT, AB, .TRUE.)
               GM = G(IM2)
               IF (I .GT. 2) THEN
                  CALL THRCOR (IM2, ISM, BET, TCOR, QTOT(IM2))
                  IF (TCOR .LE. 0.0D0) TCOR=1.0D0
                  GM = GM - RT*LOG(TCOR)
                  IF (PTEMP(IT)) WRITE (6, 614) TCOR, GM                TCA0996
               ELSE
C  for I=1,2, save it for CUS calculation
                  GEXTR(IEX,I) = GM
               END IF
C  check for absolute maximum
C              WRITE (6, 6600) ISM, SS(ISM), GM, GMAX(I)
C6600 FORMAT (' GTST: ISM,S,GM,GMAX=', I5, 1P3E13.5)
                IF (GM .GE. GMAX(I)) THEN
                  GMAX(I) = GM
                  IF (ISM .GE. NSDM1 .AND. ISM .LE. NSDM10) THEN
C  store maximum in permanent location if its in a temporary one
                     IS = NSDM + I
                     SS(IS) = SS(ISM)
                     VS(IS) = VS(ISM)
                     DS(IS) = DS(ISM)
                     XKS(IS) = XKS(ISM)
                     X2(IS) = X2(ISM)
                     Y2(IS) = Y2(ISM)
                     UXS(IS) = UXS(ISM)
                     UYS(IS) = UYS(ISM)
                     FBS(IS) = FBS(ISM)
                     QFBS(IS) = QFBS(ISM)
                     GBS(IS) = GBS(ISM)
                     XMOMS(IS) = XMOMS(ISM)
                     CAPS(IS) = CAPS(ISM)
                     ISM = IS
                  END IF
                  ISMX(I) = ISM
               END IF
C              WRITE (6, 6601) ISM, ISMX(I)
C6601          FORMAT (' GTST: ISM,ISMX=', 2I5)
  150       CONTINUE
  160       CONTINUE
            IF (I .LT. 3) THEN
               IS = ISMX(I)
C  compute kappa CVT/CAG for I=1,2
               IF (IS .LE. 0) THEN
                  IF (IOPT(4) .EQ. 0) THEN
                     E = VR
                  ELSE
                     E = VAD(1,I)
                  END IF
               ELSE IF (IS .GT. NSDM10) THEN
                  IF (IOPT(5) .EQ. 0) THEN
                     E = VP
                  ELSE
                     E = VAD(NSMAX,I)
                  END IF
               ELSE
                  E = VS(IS) + ESTR (IS, SS(IS), DS(IS), XKS(IS))
                  IF (I .EQ. 2) THEN
                     IF (LGS2(6)) THEN
C  centrifugal oscillator energy level for bend
                        CALL COBEND(NBND1,NBND2,FBS(IS),QFBS(IS),
     *                     GBS(IS), EB)
                     ELSE
C  uncoupled bending energy level
                        EB = EBEND (NBND1, FBS(IS), QFBS(IS), GBS(IS))
                        IF (NBND1 .EQ. NBND2) EB = 2.D0*EB
                        IF (NBND1 .NE. NBND2) EB = EB + EBN (NBND2,
     *                     FBS(IS), QFBS(IS), GBS(IS))
                     END IF
                     E = E + EB
                  END IF
               END IF                                                   TCA1097
               KAP(IT,2,I) = EXP(BET*(E-VMAX(I)))                       TCA1097
               VMXCVT(I) = E                                            TCA1097
C  for symmetric system, check if max in G should be reflected
               IF (LSYM .AND. NN .GT. 2) THEN
                  IS = ISMX(I)
                  IF (ABS(SS(IS)-SS(NMID)) .GT. 1.D-6) GEXTR(3,I) =
     *               GMAX(I)
               END IF
            END IF
  200    CONTINUE
C  compute and print out extra information about the transition states
         CALL EXTRAS (ISMX, IT)
         IF (PTEMP(IT)) WRITE (6, 616)                                  TCA0996
C  write out adiabatic barrier information
         DO 210 I = IC3D1,IC3D2
            IF (PTEMP(IT)) WRITE (6, 618) AH(I+1), VMXCVT(I)*CKCAL,     TCA0996
     *         VMAX(I)*CKCAL, KAP(IT,2,I)
  210    CONTINUE
C  compute rate constants
         T = -LOG(BET*TAU*TPI)
         DO 220 I = IC3D1,IC3D2
            IF (NSAD .LE. 0) THEN
               XK(IT,1,I) = 0.0D0
            ELSE
               XK(IT,1,I) = EXP(-GSP(I)/RT+T)
            END IF
            XK(IT,2,I) = EXP(-GMAX(I)/RT+T)
            XK(IT,3,I) = EXP(-GMAX(I+2)/RT+T)
            XK(IT,4,I) = XK(IT,2,I)
            IF (NEXTR(I) .EQ. 1 .OR. GEXTR(2,I) .GT. GEXTR(1,I) .OR.
     *         GEXTR(2,I) .GT. GEXTR(3,I)) GO TO 220
            IF (GMAX(I) .EQ. GEXTR(3,I)) THEN
               I2 = 1
            ELSE
               I2 = 3
            END IF
            R = 1.0D0 + EXP((GEXTR(I2,I)-GMAX(I))/RT) - EXP((GEXTR(2,I)-
     *         GMAX(I))/RT)
            XK(IT,4,I) = XK(IT,2,I)/R
  220    CONTINUE
C
         IF (PTEMP(IT)) WRITE (6, 622)                                  TCA0996
         IF (LGS(9) .AND. PTEMP(IT)) WRITE (6, 624) (XK(IT,J,1),J=1,4)  TCA0996
         IF (LGS(10) .AND. PTEMP(IT)) WRITE (6, 626) (XK(IT,J,2),J=1,4) TCA1296
  300 CONTINUE
      RETURN
  600 FORMAT (/, 1X, 22('*'), 2X, 'CVT, ICVT, and CUS Calculations ',
     *               22('*'))
  602 FORMAT (/, 1X, T5, 'Free energy curves, G(s), in units of ',
     *                   'kcal (standard state is ',
     *        /, 1X, T5, 'molec/cm (collinear), molec/cc (3D))',
     *        /, 1X, T5, 'All rotational partition functions ',
     *                   'include the symmetry factors',
     *        /, 1X, T5, A9, ' vibrational partition functions ',
     *                   'used for stretch')
  604 FORMAT (1X, T5, 'Harmonic-Quartic bend potential used')
  605 FORMAT (1X, T5, 'Harmonic bend potential used')
  606 FORMAT (//, 1X, 10('*'), 1X,'Temperature = ', F8.2, 
     *            5X, 'Beta = ', 1PE13.5, ' Hartree**-1',
     *        //, 1X, T5,'s',T15,'ln',T24,'Rot', T28,
     *   'Quantum', T36, 'L**2/(4PI)', T47, 'delta PHI', T60, 'ln',
     *   T68, 'Quartic', T78, 'ln', T86, A9, T97, 'G(s)', T111,
     *   'G(s)',/,T14,'QROT',T24,'Sig',T29,'Corr.',T48,'Degrees',
     *   T57, 'Q(BND)**2', T69, 'Corr.', T77, 'Q(STR)', T88, 'Corr.',
     *   T95, 'Collinear', T112, '3D'/)
  607 FORMAT (1X, T5, 'The collinear asymptotic reactant value ',
     *                'is not used in the search',
     *     /, 1X, T5, 'for the maximum in G.')
  608 FORMAT (/, 1X, T5, 'For symmetric reactions the results are ',
     *                   'reflected across the symmetric',
     *        /, 1X, T5, 'stretch line.')
  609 FORMAT (1X, T5, 'The collinear asymptotic product value is ',
     *                'not used in the search for the ',
     *     /, 1X, T5, 'maximum in G.')
  610 FORMAT (/, 1X, T5, 'Extrema in G(s,T)', 2X, A15, /)
  612 FORMAT (/, 1X, T5, 'Extrema in ICVT curves', 2X, A15,
     *        /, 1X, T5, 'The free energies do not contain ',
     *                   'the threshold correction', /)
  614 FORMAT (1X, T5, 'Threshold correction = ', 1PE13.5, 
     *                ', G(corrected) = ',1PE13.5)
  616 FORMAT (/, 1X, 20('*'),1X,'Adiabatic barrier heights (kcal)',
     *        /,26X,'CVT',11X,'AG',5X,'Kappa CVT/CAG')
  618 FORMAT (1X, A15, 5X, 1P,3E13.5)  
  620 FORMAT (/, 1X, T5, 'Maxima in the polarization potential ',
     *                   'are at ', 3F12.6, 
     *        /, 1X, T5, 'for E = KT, E = 5*KT/3, and E = 2KT, ',
     *                   'respectively.' )
  622 FORMAT (/, 1X, 20('*'), 1X, 'Rate constants',/,22X,'TST',11X,     TCA0197
     *   'CVT',11X, 'ICVT', 10X, 'CUS', 10X, 'Units')                   TCA0197
  624 FORMAT (' Collinear', 6X, 1P,4E14.5, 5X, 'cm/molec-sec')          TCA0197
  626 FORMAT (' 3D', 13X, 1P,4E14.5, 2X, 'cc/molec-sec')     
  630 FORMAT (//,1X,10('*'),1X,'Temperature = ',F8.2,' free energy ',   TCA0996
     *        'data not printed for this temperature',/)                TCA0996
 6000 FORMAT (/,2X,T5,'Iteration to find GMAX not converged, DELS=',
     *   1PE13.5, ', last 3 G=', 1P,3E13.5)                 
      END
