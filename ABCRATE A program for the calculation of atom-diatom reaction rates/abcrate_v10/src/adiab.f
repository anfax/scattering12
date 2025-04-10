!!!*************************************************************
! 文件/File: adiab.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: adiab.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:40
!*************************************************************

      SUBROUTINE ADIAB
C
C     ADIAB  - compute adiabatic potential curves
C
C     Modified 3/18/91 to include centrifugal oscillator bend energies
C
C  Called by:
C     VTSTM  - main program
C
C  Calls:
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C     DERS   - derivatives of morse turning point and zero pt. energy
C              w.r.t. s
C     EBEND  - compute bending energy levels
C     ESTR   - compute stretch vibrational energy levels
C     INTERP - interpolate r.p. info from grid
C     QUADFT - quadratic fit of three points
C     SAGCMP - compute info needed for effective mass terms
C     STORE  - store reaction path info
C     TABL21 - print out table of GTS info
C     TITLES - print out 2-line title
C     WKBSET - set up grid of WKB energy levels
C     WKBWRT - write out WKB energies, etc. to unit 14
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*5 LEVTYP(3)
      CHARACTER*2 C1D3D(2)
      LOGICAL LCONV, LPRNT
      DIMENSION V(2), VSV(4,2), B(3), VASP(2), VMAX2(2), ISMAX2(2)
C  adiabatic barriers are stored in VAD(I,J), I-index over s values,
C    J=1 is collinear, J=2 is 3d
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /ADIAB1/ VAD(NSDM,2), VMAX(2), VR, VP, SMAX(2), ISMAX(2)
      COMMON /ADIAB2/ SAD1, SAD2, NPAD
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /EBND2/  ALF, ALFBAR
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      LOGICAL LGS2(10)
      COMMON /LOGIC2/ LGS2
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /MORAB/  DAB, XKAB, AMAB, VDELTA
      COMMON /MORBC/  DBC, XKBC, AMBC
      COMMON /OPTION/ IOPT(20)
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     *XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     *NSAD, NSADMX(2)
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      COMMON /STATE2/ DGBND, LSBEND, NBND1, NBND2
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL1096
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
      SAVE LEVTYP, C1D3D                                                TCA1097
      DATA LEVTYP /'Morse', ' WKB ', 'Harm.'/
      DATA C1D3D /'1D', '3D'/
      NSLAST = NSMAX
      IF (LSYM) NSLAST = NMID
      N = NSADMX(2)
      NSHLF = NSSP(N)
C
      IF (LGS(7)) THEN
C  compute WKB energy levels for stretch
         ILTYP = 2
         IF (.NOT.LGS2(5)) THEN
            DO 10 IS = 1,NSLAST
               S = SS(IS)
               D = DS(IS)
               XKM = XKS(IS)
C  give Morse eigenvalue as initial guess for state NSTATE
               E = D*TNP1*(2.D0-TNP1/XKM)/XKM
               CALL WKBSET (IS, S, NSTATE, E, UL, UG, PER)
C  if state-selected and not for the ground state compute ground state
C     WKB eigenvalue
               IF (LSTATE .NE. 0  .AND. NSTATE .NE. 0) THEN
                  E = D*(2.D0-1.D0/XKM)/XKM
                  CALL WKBSET (IS, S, 0, E, UL, UG, PER)
               END IF
   10       CONTINUE
            IF (LGS2(4)) CALL WKBWRT
         END IF
      ELSE
         IF (LGS(3)) THEN
            ILTYP = 1
         ELSE
            ILTYP = 3
         END  IF
      END IF
C  write out information about generalized transition at saddle point
C     to lfn 21
      CALL TABL21 (0.D0, NSHLF, 0, 1)
      ISMAX(1) = 0
      ISMAX(2) = 0
      SMAX(1) = -1.D30
      SMAX(2) = -1.D30
C  reactant adiabatic energy (assumes bends and rotations in ground
C     state)
      VR = ESTR (-1, -1000.D0, DBC, XKBC)
      WRITE (6, 600)
      CALL TITLES (1, 6, 1)
      WRITE (6, 601) NSTATE, NBND1, NBND2, LEVTYP(ILTYP)
      IF (IOPT(4) .EQ. 0) THEN
C  include reactant values in search for maximum in VA
         VMAX(1) = VR
         VMAX(2) = VR
         VSV(4,1) = VR
         VSV(4,2) = VR
      ELSE
C  don't include reactant values in search for maximum in VA
         WRITE (6, 602)
         VMAX(1) = -1.D30
         VMAX(2) = -1.D30
         VSV(4,1) = -1.D30
         VSV(4,2) = -1.D30
      END IF
      T1 = VR*CKCAL
      WRITE (6, 603) T1, T1
      ISINT = 1
C  loop over grid of stored values, compute and store VA, locate max
C     on this grid
      IPRNT = NPAD
      ISP = 0
      ISAD = 0
      DO 30 IS = 1,NSLAST
C  check if information should be printed, if so set LPRNT = .true.
C     print forced at all saddle points
         IF (IS .GE. ISP) THEN
            ISAD = ISAD + 1
            IF (ISAD .LE. NSAD) THEN
               ISP = MIN(NSLAST, NSSP(ISAD))
            ELSE
               ISP = NSLAST
            END IF
            LPRNT = .TRUE.
         ELSE
            IPRNT = IPRNT + 1
            LPRNT = IPRNT .GE. NPAD
         END IF
         IF (LPRNT) IPRNT = 0
         S = SS(IS)
         XKM = XKS(IS)
         D = DS(IS)
C  collinear VA (note ESTR will return WKB energy if LGS(7) is true)
         V(1) = VS(IS) + ESTR (IS,S,D,XKM)
C  bending potential, allowed to be state selected to be consistent
C     with RESON program at CDC.
         IF (LGS2(6)) THEN
C  centrifugal oscillator bending energy level
            CALL COBEND(NBND1,NBND2,FBS(IS),QFBS(IS),GBS(IS), EB)
         ELSE
C  uncoupled bending energy level
            EB = EBEND (NBND1, FBS(IS), QFBS(IS), GBS(IS))
            IF (NBND1 .EQ. NBND2) THEN
               EB = 2.D0*EB
            ELSE
               EB = EB + EBN (NBND2, FBS(IS), QFBS(IS), GBS(IS))
            END IF
         END IF
         V(2) = V(1) + EB
         VAD(IS,1) = V(1)
         VAD(IS,2) = V(2)
         IF (LPRNT) THEN
C  mu/mubar factors for MCPSAG, SCSAG and PASAG computed for output
C     only
            CP = CAPS(IS)
            CALL SAGCMP (S, D, XKM, CP, AM, OM, XLM, UM)
            CALL DERS (ISINT, S, DUDSM, DEDS)
            CPUM = CP*UM
            D2 = DUDSM*DUDSM
            A1 = 1.D0 - CPUM
            IF (A1 .LT. 0.0D0) A1 = 0.0D0
            A1 = A1*A1 + D2
            T = -CPUM*(2.D0+CPUM) + D2
            A2 = 1.0D0
            IF (T .LT. 0.0D0) A2 = EXP(T)
            T = 1.0D0 - XLM
            A3 = 0.D0                                                   GCL1092
            IF (T .GT. 0.0D0) THEN
               A3 = SQRT(T)*T
               IF (LGS(3)) A3 = T /(1.D0+3.D0*CP*TNP1/(AM*XKM))
            END IF
            T1 = V(1)*CKCAL
            T2 = V(2)*CKCAL
            IF (NSAD .GT. 0 .AND. IS .EQ. NSHLF) WRITE (6,604)
C  alf and alfbar are effective harmonic frequencies needed to
C     reproduce the harmonic and harm.-quartic bend energy level for
C     n=0.  These are computed in EBEND
            WRITE (6, 606) SS(IS), T1, T2, CP, UM, A2, ALF, ALFBAR
            IF (NSAD .GT. 0 .AND. IS .EQ. NSHLF) WRITE (6,604)
C
C           WRITE (80, 8000) X2(IS), Y2(IS)
C           WRITE (81, 8000) X2(IS)+UXS(IS)*UM,Y2(IS)+UYS(IS)*UM
C8000 FORMAT (1X, 1P2E20.10)
C
         END IF
C  check if the  maximum
         DO 20 I = 1,2
            IF (IS .EQ. NSHLF) VASP(I) = V(I)
            IF (V(I) .GE. VMAX(I)) THEN
               ISMAX(I) =IS
               VMAX(I) = V(I)
C  VSV holds values around maximum for use in quadratic fit below.
               VSV(1,I) = VSV(4,I)
            END IF
            IF (IS-ISMAX(I) .EQ. 1) VSV(3,I) = V(I)
            VSV(4,I) = V(I)
   20    CONTINUE
         IF (LSYM) THEN
C  if symmetry option on VA reflected across symmetry point
            I2 = 2*NMID-IS
            IF (I2 .LE. NSMAX .AND. I2 .NE. NMID) THEN
               VAD(I2,1) = VAD(IS,1)
               VAD(I2,2) = VAD(IS,2)
            END IF
         END IF
   30 CONTINUE
C  loop over grid points done, now set product value of VA
      IF (LSYM) THEN
C  symmetric system, product value equals reactant value
         VP = VR
         WRITE (6, 608)
      ELSE
C  not symmetric, compute product value
         VP = VDELTA + ESTR (-2, 1000.D0, DAB, XKAB)
         T1 = VP*CKCAL
         WRITE (6, 610) T1, T1
         IF (IOPT(5) .EQ. 0) THEN
C  check if product value is maximum
            DO 40 I = 1,2
               IF (VMAX(I) .GT. VP) GO TO 40
               ISMAX(I) = 0
               SMAX(I) = 1.D30
               VMAX(I) = VP
   40       CONTINUE
         ELSE
C  don't check if product value is maximum
            WRITE (6, 611)
         END IF
      END IF
C  find max. of adiabatic curve (ISMAX contains locations of max. on
C     the grid).  Loop over collinear and 3D.
      DO 70 I = 1,2
         IS = ISMAX(I)
         IF (LSYM .AND. IS .EQ. NSHLF) THEN
C  if symmetric system and maximum at saddle point don't do quadratic
C     fit
            S = 0.0D0
         ELSE IF (IS .LE. 0) THEN
C  if maximum at reactants or products (IS=0)
            S = SMAX(I)
         ELSE
            S = SS(IS)
C  if maximum at either end of the grid skip quadratic fit
            IF (IS .NE. 1 .AND. IS .NE. NSLAST) THEN
C  Quadratic fit is performed on the maximum and points on either side
C     of it.  The maximum of this quadratic fit is found and the
C     potential parameters and adiabatic potential are found for it and
C     a point on either side.  These three points are then used for the
C     quadratic fit.  The step size for the points on either sides is
C     halfed if the predict maximum lies within the previous
C     three-point grid, otherwise it is unchanged.  The procedure is
C     terminated if the predicted and actual maximum values are within
C     1.E-6 (relative error), the previous and newly computed VA are
C     within 1.E-6, the step size is below 1.E-9, or 10 iterations have
C     been made.
C
C  ISI - location at the first point in three-point grid in the large
C     grid.  (Used for interpolation)
C  ISST - location in large grid for storeage of the newly computed
C     potential information.
               ISI = IS - 1
               ISST = NSDM4
               ISSTD = 3
               DELS = 0.5D0 * (SS(IS) - SS(ISI))
               VSV(2,I) = VMAX(I)
               JMAX = 2
               SMX = S
               II = 0
C  Loop over iterations
   50          CONTINUE
                  II = II + 1
                  VM = 0.D0                                             GCL1092
                  S = SMX
                  IF (JMAX .EQ. 2) THEN
C  quadratic fit of VAD around maximum
                     CALL QUADFT (SS(ISI), VSV(1,I), B)
                     IF (B(3) .LT. 0.D0) THEN
                        S = -.5D0*B(2)/B(3)
                        IF (S .GT. SS(ISI) .AND. S .LT. SS(ISI+2))
     *                     THEN
                           VM = B(1)+S*(B(2)+S*B(3))
                        ELSE
                           S = SMX
                        END IF
                     END IF
                  END IF
                  S = S - 2.0D0*DELS
                  JMAX = 0
                  LCONV = .TRUE.
C  Compute potential parameters and VA on three new points and store
C     them in the grid.
                  DO 60 J = 1,3
                     S = S + DELS
C  First interpolate X,Y to the point S from the large grid and compute
C     potential parameters at that X,Y location
                     CALL INTERP (ISI)
C  Store the potential information
                     IS = ISST + J
                     CALL STORE (IS)
C  Compute adiabatic potential
                     VSV(J,I) = VS(IS) + ESTR(IS, S, DS(IS), XKS(IS))
                     IF (I .EQ. 2) THEN
                        IF (LGS2(6)) THEN
                           CALL COBEND (NBND1,NBND2,FBS(IS),QFBS(IS),
     *                        GBS(IS),EB)
                        ELSE
                           EB = EBEND(NBND1, FBS(IS), QFBS(IS), GBS(IS))
                           IF (NBND1 .EQ. NBND2) THEN
                              EB = 2.D0*EB
                           ELSE
                              EB = EB + EBN(NBND2, FBS(IS), QFBS(IS),
     *                         GBS(IS))
                           END IF
                        END IF
                        VSV(J,I) = VSV(J,I) + EB
                     END IF
C  See if this point is the maximum
                     IF (VSV(J,I) .GE. VMAX(I)) THEN
                        JMAX = J
                        SMX = S
                        VMAX(I) = VSV(J,I)
                     END IF
                     IF (J .GT. 1) LCONV = LCONV .AND.
     *                  ABS(1.0D0-VAOLD/VSV(J,I)) .LT. 1.D-6
                     VAOLD = VSV(J,I)
   60             CONTINUE
                  IF (.NOT.LCONV) THEN
C  If maximum is on three-point gird reset ISI and ISST so that the
C     three new points are used for interpolation and the next three
C     points are stored at a new location in the grid.
                     IF (JMAX .NE. 0) THEN
                        ISI = ISST + 1
                        ISST = ISST + ISSTD
                        ISSTD = -ISSTD
                        LCONV = ABS(1.0D0-VM/VMAX(I)) .LT. 1.D-6
                     END IF
C  Halve step size if maximum was at middle or none of the new values
C     beat the old maximum.
                     IF (JMAX .EQ. 0 .OR. JMAX .EQ. 2) DELS = 0.5D0*DELS
                     LCONV = LCONV .OR. DELS .LT. 1.D-9
                  END IF
               IF (II .LT. 10 .AND. .NOT.LCONV) GO TO 50
C  End of loop over iterations in quadratic fits
               S = SMX
               IF (ABS(S) .LE. 1.D-4) THEN
C  Maximum so close to saddle point that saddle point is assumed.
               ISMAX(I) = NSHLF
               S = 0.0D0
            ELSE
               IF (.NOT.LCONV) WRITE (6, 6000) I, DELS,
     *            (VSV(J,I)*CKCAL,J=1,3)
C  Interpolate and store potential information at the maximum.
                  CALL INTERP (ISI)
                  IS = NSDM + I
                  CALL STORE (IS)
                  ISMAX(I) = IS
               END IF
            END IF
         END IF
         SMAX(I) = S
C  write out information about generalized transition at adiabatic
C     maximum to lfn 21
         CALL TABL21 (0.D0, IS, 1, I)
   70 CONTINUE
C  end of loop over 1d and 3d
      DO 80 I = 1,2
         VSV(1,I) = VMAX(I)*CKCAL
         VSV(2,I) = VASP(I)*CKCAL
         VSV(3,I) = (VASP(I)-VMAX(I))*CKCAL
   80 CONTINUE
      WRITE (6, 612) (SMAX(I), (VSV(J,I), J=1,3), I=1,2)
C
      IF (SAD1 .GE. SAD2) RETURN
C  Find another local maximum between SAD1 and SAD2
      IS1 = 1
C  Find limits on grid index
      DO 90 IS = 1,NSMAX
         IF (SS(IS) .GT. SAD1) GO TO 100
         IS1 = IS
   90 CONTINUE
  100 CONTINUE
      IF (IS .EQ. NSMAX) RETURN
      DO 110 IS = IS1,NSMAX
         IS2 = IS
         IF (SS(IS) .GT. SAD2) GO TO 120
  110 CONTINUE
  120 CONTINUE
      IF (IS2 .LE. IS1) RETURN
      VMAX2(1) = -1.D10
      VMAX2(2) = -1.D10
C  Loop over grid points between SAD1 and SAD2
      DO 140 IS = IS1,IS2
         DO 130 I = 1,2
            IF (VAD(IS,I) .LT. VMAX2(I)) GO TO 130
            ISMAX2(I) = IS
            VMAX2(I) = VAD(IS,I)
  130    CONTINUE
  140 CONTINUE
C  Simple quadratic fits
      DO 150 I = 1,2
         IF (ISMAX2(I) .EQ. IS1 .OR. ISMAX2(I) .EQ. IS2) GO TO 150
         ISI = ISMAX2(I) - 1
         ISI = MAX (1, ISI)
         ISI = MIN (ISI, NSMAX-2)
         CALL QUADFT (SS(ISI), VAD(ISI,I), B)
         S = -.5D0*B(2)/B(3)
         VM = B(1) + S*(B(2) + S*B(3))
         WRITE (6, 614) C1D3D(I), SAD1, SAD2, S, VM*CKCAL
  150 CONTINUE
      RETURN
C
  600 FORMAT (/,1X, 25('*'),2X,'Adiabatic barriers (kcal)',1X,25('*'))
  601 FORMAT (/,' NSTR=',I5,', NBEND-',2I5,/,T64,'__',T76,'Bend ALPHA ',TCA0197
     *   '(unitless)',/,17X,'Adiab. curves',23X,A5,5X,'MU/MU',7X,       TCA0197
     *   24('-'),/,5X,'s',7X,'Collinear',6X,'3D',9X,'Curvature',3X,     TCA0197
     *   'Turn. pt.',3X,'CD-SCSAG',7X,'Harm.',8X,'H.-Q.')               TCA0197
  602 FORMAT (/,1X,T5,'The reactant values are not used in the search ',
     *                'for the maximum',
     *        /,1X,T5,'in the vibrational adiabatic barriers')
  603 FORMAT (/,1X,T3,'Reactants', 1P,2E13.5,/)                         TCA0197
  604 FORMAT (' ')
  606 FORMAT (1X, F8.4, 1X, 1P,3E13.5, 0PF10.5, 1PE13.5, 1X, 2E13.5)
  608 FORMAT (/,1X,T5,'Symmetric system, parameters on the product',
     *                ' side of the symmetric stretch',
     *        /,1X,T5,'line are found by reflection')
  610 FORMAT (/,1X,T5,'Products', 1X, 1P,2E13.5) 
  611 FORMAT (/,1X,T5,'The product values are not used in the search ',
     *                'for the maximum',
     *        /,1X,T5,'in the vibrational adiabatic barriers')
  612 FORMAT (/,1X,T5, 'Adiabatic barrier heights (kcal)',
     *        /,1X,T5, 'Collinear', T16,'s(max) = ', 1PG12.5,           TCA1097
     *            T40, 'Va(max)', T55, '=', T57, E15.7,                 TCA1097
     *       /,1X,T40, 'Va(saddle pt.)', T55,'=', T57, 1PE15.7, 
     *       /,1X,T40, 'Difference', T55, '=', T57, 1PE15.7,
     *        /,1X,T5, '3D', T16, 's(max) = ', 1PG12.5,                 TCA1097
     *            T40, 'Va(max)', T55, '=', T57, E15.7,                 TCA1097
     *       /,1X,T40, 'Va(saddle pt.)', T55, '=', T57, 1PE15.7, 
     *       /,1X,T40, 'Difference', T55, '=', T57, 1PE15.7)
  614 FORMAT (/,1X,T5, 'Local maximum in ', A2, ' curve between s = ', 
     *        F15.8, ' and ', F15.8, 
     *        /,1X,T5, 'occurs at s = ', F15.8, ', VaD = ', 1PE13.5,
     *                 ' kcal')
 6000 FORMAT(/,1X,T5,'Iteration to Va(max) not converged for I = ', I2, 
     *       ', DELS = ', 1PE13.5, 
     *       /,1X,T5,'The last 3 values of Va = ', 1P, 3E13.5)
      END
