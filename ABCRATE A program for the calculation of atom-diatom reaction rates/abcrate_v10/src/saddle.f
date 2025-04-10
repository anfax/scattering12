!!!*************************************************************
! 文件/File: saddle.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: saddle.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE SADDLE
C
C     SADDLE - find saddle point(s) and do normal mode analysis
C
C     Modified 3/18/91 to include centrifugal oscillator bend energies
C
C  Called by:
C     DATAIN - read in data and set up constants
C
C  Calls:
C     BEND   - compute bending potential parameters
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C     D2DX2  - compute matrix of second derivatives in Jacobi coor.
C     PEF    - evaluate potential
C     ROOT2D - 2d root search
C     RSP2   - diagonalize 2x2 matrix
C     WKBINT - extrapolate WKB energy levels from grid
C     WKBSET - set up grid of WKB energy levels
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL DERIV
      LOGICAL LP
      DIMENSION F(3), ROOT(2), R1S(11), SCR(9), SCR2(9), SCR3(9),
     *SCR4(9), VAD(2), VECT(2,2), VMX(2), ISCR(11)
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /BENDTS/ FBTS(NSADDM), QFBTS(NSADDM), GBTS(NSADDM),
     *XMOMTS(NSADDM)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      LOGICAL LMAX
      COMMON /EBND1/   LMAX
      COMMON /EBND2/  ALF, ALFBAR
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      LOGICAL LGS2(10)
      COMMON /LOGIC2/ LGS2
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /MORBC/  DBC, XKBC, AMBC
      COMMON /MORSTS/ DTS(NSADDM), XKTS(NSADDM), AMTS(NSADDM),
     *OMTS(NSADDM), OMIMG(NSADDM)
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     *XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     *NSAD, NSADMX(2)
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      PARAMETER (NLEVEL=9)
      DIMENSION NN(NLEVEL),NOM(NLEVEL),DEGEN(NLEVEL)
      EQUIVALENCE (LP, LGS(1))
      SAVE NN, NOM, DEGEN                                               TCA1097
      DATA NN   /0,  0,  0,  1,  0,  1,  0,  1,  2/
      DATA NOM  /0,  1,  2,  0,  3,  1,  4,  2,  0/
      DATA DEGEN/1.0D0,2.0D0,2.0D0,1.0D0,2.0D0,2.0D0,2.0D0,2.0D0,1.0D0/ GCL0992
C
C  Loop over saddle points. SADDLE is only called if NSAD > 0.
      VMX(1) = -1.D30
      VMX(2) = -1.D30
      DO 100 ISP = 1,NSAD
C  R1SP, R2SP are the initial guess of the saddle point geometry
C  R1, R2 are the A-B, B-C distances, respectively.
         IF(LP) WRITE (6, 600) ISP,R1SP(ISP), R2SP(ISP)                 TCA0497
         X = R1SP(ISP) + CM1*R2SP(ISP)
         Y = CM2*R2SP(ISP)
         NMAX = 25
C  root search for saddle point (solve for zero derivatives)
         CALL ROOT2D(DERIV, X, Y, 1.D-4, 1.D-10, 1.D-8, NMAX)
         IF (NMAX .EQ. 0) THEN
C  root not found, print error message and grid of V(R1,R2)
            WRITE (6, 6000)
            R2 = R2SP(ISP) + 0.6D0
            DO 20 I = 1,11
               R2 = R2 - 0.1D0
               R1 = R1SP(ISP) - 0.6D0
               DO 10 J = 1,11
                  R1 = R1 + 0.1D0
                  R3 = R1 + R2
                  CALL PEF (R1, R2, R3, V, D1, D2, D3, 0)               GCL0893
                  VV = V
                  ISCR(J) = INT(1000.D0*CKCAL*VV + 0.5D0)
   10          CONTINUE
               WRITE (6, 6001) R2, ISCR
   20       CONTINUE
            R1 = R1SP(ISP) - 0.6D0
            DO 30 I = 1,11
               R1 = R1 + 0.1D0
               R1S(I) = R1
   30       CONTINUE
            WRITE (6, 6002) R1S
            STOP 'SADDLE 1'
         END IF
C  evaluate potential at saddle point
         XSP(ISP) = X
         YSP(ISP) = Y
         R2 = YSP(ISP)/CM2
         R1 = XSP(ISP) - CM1*R2
         R3 = R1 + R2
         CALL PEF (R1, R2, R3, V, D1, D2, D3, 1)                        GCL0893
         VSP(ISP) = V
         IF (LP) THEN
            T1 = VSP(ISP)*EAU
            T2 = VSP(ISP)*CKCAL
            DX = D1 + D3
            DY = (D2+D3 - CM1*DX)/CM2
            WRITE (6, 602) R1, R2, VSP(ISP), T1, T2, D1, D2, D3, DX, DY
         END IF
         R1SP(ISP) = R1
         R2SP(ISP) = R2
C  compute second derivative matrix
         CALL D2DX2 (XSP(ISP), YSP(ISP), F, ERR)
         IF (LP) WRITE (6, 604) F, ERR
C  diagonalize second derivative matrix
         CALL RSP2 (F, ROOT, 1, VECT)
C  check for saddle point
         IF ((ROOT(1)*ROOT(2)) .GE. 0.D0) THEN                          GCL1092
            WRITE (6, 6003)
            STOP 'SADDLE 2'
         END IF
C  normal modes
         DELRBC = VECT(2,1)/CM2
         DELRAB = VECT(1,1) - CM1*DELRBC
         IF (DELRBC .LT. DELRAB) THEN
            SVECT(1,ISP) = -VECT(1,1)
            SVECT(2,ISP) = -VECT(2,1)
         ELSE
            SVECT(1,ISP) = VECT(1,1)
            SVECT(2,ISP) = VECT(2,1)
         END IF
         IF (SVECT(1,ISP)*VECT(2,2) .GT. 0.0D0) THEN
            UVECT(1,ISP) = -VECT(1,2)
            UVECT(2,ISP) = -VECT(2,2)
         ELSE
            UVECT(1,ISP) = VECT(1,2)
            UVECT(2,ISP) = VECT(2,2)
         END IF
         IF (LP) WRITE (6, 606) ROOT, (SVECT(I,ISP), UVECT(I,ISP),
     *      I=1,2)
C  transition state morse curve parameters
         OMTS(ISP) = SQRT(ROOT(2)/XMU)
         OMIMG(ISP) = SQRT(-ROOT(1)/XMU)
         DTS(ISP) = DBC - VSP(ISP)
         XKTS(ISP) = 4.D0*DTS(ISP)/OMTS(ISP)
         AMTS(ISP) = 2.D0*SQRT(2.D0*XMU*DTS(ISP))/XKTS(ISP)
         VAD(1) = VSP(ISP) +
     *      DTS(ISP)*OMTS(ISP)*0.5D0*(1.D0-OMTS(ISP)/XKTS(ISP))
         IF (LP) THEN
C   write out harmonic frequencies
            WRITE (6, 608) OMTS(ISP), OMTS(ISP)*EAU, OMTS(ISP)*CKCAL,
     *         OMTS(ISP)*CCM, OMIMG(ISP), OMIMG(ISP)*EAU,
     *         OMIMG(ISP)*CKCAL, OMIMG(ISP)*CCM
            XE = 1.D0/XKTS(ISP)
            T1 = DTS(ISP)*EAU
            T2 = DTS(ISP)*CKCAL
            WRITE (6, 610) DTS(ISP), T1, T2, AMTS(ISP), XKTS(ISP), XE
         END IF
C   WKB eigenvalues
         IF (LGS(7)) THEN
            IS = -2-ISP
            IF (.NOT.LGS2(5)) THEN
               E = 0.5D0*OMTS(ISP)*(1.D0-0.5D0*XE)
               S = 0.D0                                                 GCL1092
               CALL WKBSET (IS, S, 0, E, UL, UG, PER)
               VAD(1) = VSP(ISP) + E
            ELSE
               CALL WKBINT (IS, S, 0, E, UL, UG, PER)
            END IF
            IF (LP) WRITE (6, 612) E, E*EAU, E*CKCAL, E*CCM
            IF (LSTATE .NE. 0 .AND. NSTATE .NE. 0) THEN
               IF (.NOT.LGS2(5)) THEN
                  E = 0.5D0*TNP1*OMTS(ISP)*(1.D0-0.5D0*TNP1*XE)
                  CALL WKBSET (IS, S, NSTATE, E, UL, UG, PER)
               ELSE
                  CALL WKBINT (IS, S, NSTATE, E, UL, UG, PER)
               END IF
               IF (LP) WRITE (6, 614) NSTATE, E, E*EAU, E*CKCAL, E*CCM
            END IF
         END IF
         IF (LP) THEN
            WRITE (6, 616)
C   energy levels of saddle point stretch
            NMAX = INT(0.5D0*(XKTS(ISP)-1.0D0)) + 1
            TIXK = 2.0D0/XKTS(ISP)
            T = -0.5D0*TIXK
            DO 40 N = 1,NMAX
               NQM = N-1
               T = T + TIXK
               E = DTS(ISP)*T*(2.D0-T)
               WRITE (6, 618) NQM, E, E*EAU, E*CKCAL, E*CCM
   40       CONTINUE
         END IF
         CALL BEND (XSP(ISP), YSP(ISP), FB, QFB, GB, XMOM)
         FBTS(ISP) = FB
         QFBTS(ISP) = QFB
         GBTS(ISP) = GB
         XMOMTS(ISP) = XMOM
C
         IF (LGS2(6)) THEN
C  Centrifugal oscillator energy levels for bend
            CALL COBEND(0,0,FB,QFB,GB,EB)
            VAD(2) = VAD(1) + EB
         ELSE
C  uncoupled bending energy level
            EB = EBEND(0,FB,QFB,GB)
            VAD(2) = VAD(1) + 2.0D0*EB                                  GCL1092
         END IF
         IF (LP) THEN
            IF (FB .LE. 0.D0) THEN
               IF (.NOT.LGS(4)) WRITE (6, 620)
               T1 = 0.D0                                                GCL1092
            ELSE
               T1 = SQRT(FB*GB)*CCM
            END IF
            WRITE (6, 622) FB*CKCAL, QFB*CKCAL, T1, GB, XMOM
            W = T1*CKCAL/CCM
            IF (LGS2(6)) THEN
C  Centrifugal oscillator energy levels for bend
               I = 1
               T =  W
               SCR(1) = T
               SCR2(1) = EB*CKCAL
               DO 45 N = 1,4
                  T = (N+1)*W
                  K = N + 2
46                CONTINUE
                     K = K - 2
                     NV = (N-K)/2
                     CALL COBEND(NV,K,FB,QFB,GB,EB)
                     IF (.NOT.LMAX) THEN
                        I = I + 1
                        SCR(I) = T
                        SCR2(I) = EB*CKCAL
                     END IF
                  IF (K.GT.1) GO TO 46
45             CONTINUE
               N = I
               IF (FB .LE. 0.D0) WRITE (6, 624)
               WRITE (6, 626) (SCR(I), I=1,N)
               WRITE (6, 627) (SCR2(I), I=1,N)
C               
            ELSE
C  uncoupled bending energy levels
               T = 0.5D0*W                                              GCL1092
               SCR(1) = T
               SCR2(1) = EB*CKCAL
               SCR3(1) = ALF
               SCR4(1) = ALFBAR
               N = 0
               DO 50 I = 2,9
                  T = T + W
                  SCR(I) = T
                  N = N + 1
                  SCR2(I) = EBN(N, FB, QFB, GB)*CKCAL
                  IF (LMAX) THEN
                     N = N - 1
                     GO TO 60
                  END IF
                  SCR3(I) = ALF
                  SCR4(I) = ALFBAR
50             CONTINUE
60             CONTINUE
               N =N +1
               IF (FB .LE. 0.D0) WRITE (6, 624)
               WRITE (6, 628) (SCR(I), I=1,N)
               WRITE (6, 627) (SCR2(I), I=1,N)
               IF (FB .GT. 0.D0) THEN
                  IF (LGS(4)) THEN                                      GCL1296
                      IF (LGS(8)) THEN
                          WRITE (6, 630)
                      ELSE
                          WRITE (6, 632) (SCR3(I), I=1,N)
                          WRITE (6, 628) (SCR4(I), I=1,N)
                      END IF
                  END IF                                                GCL1296
               ELSE
                  WRITE (6, 634)
               END IF
            END IF
         END IF
         DO 70 IC3D = 1,2
            IF (VAD(IC3D) .LT. VMX(IC3D)) GO TO 70
            IF (ABS(VAD(IC3D)-VMX(IC3D)) .LT. 1.D-6) GO TO 70
            NSADMX(IC3D) = ISP
            VMX(IC3D) = VAD(IC3D)
   70    CONTINUE
  100 CONTINUE
      RETURN
  600 FORMAT (/,1X,26('*'),2X,'Saddle point properties',1X,26('*'),//,
     *   1X,'Saddle point ',I1,/,1X,T5,'Geometry',/,1X,T15,'Initial ',  TCA0197
     *   'guess',T30,'RAB=',F12.8,T50,'RBC=', F12.8, T70, 'bohr')
  602 FORMAT (1X,T15,'Exact',T30,'RAB=', F12.8,T50,'RBC=',F12.8,T70,
     *   'bohr', //,1X,T5,'Potential energy, VSP=',T30,1PE17.8, 
     *   T50,'Hartree',T60,1PE17.8,T78,'eV',/,1X,T30,1PE17.8,T50,
     *   'kcal',//,1X,T5,'Derivatives (a.u.)',T25,'DR1=',T32,1PE17.8,
     *   T55,'DR2=',T62,1PE17.8,/,1X,T25,'DR3=',T32,1PE17.8,
     *   /,1X,T25,'DX=',T32,1PE17.8,T55,'DY=',T62,1PE17.8)
  604 FORMAT (/,1X,T5,'Force constants (FXX,FXY,FYY) =',T40,1PE17.8, 
     *   T60,1PE17.8,/,1X,T40,1PE17.8,T60,'Hartree/bohr**2',
     *   /,1X,T5,'Average relative error is', E13.5,
     *   //,1X,T5,'Coordinates (X,Y) defined by', T40, 'X = RAB +',
     *   ' (C/(B+C))*RBC',/,1X,T40,'Y = SQRT(M/MU)*RBC')
  606 FORMAT (/,1X,T5,'Eigenvalues of F',T30,1PE20.10,T50,1PE20.10,
     *   /,1X,T5,'Eigenvectors of F',(T30,1PE20.10,T50,1PE20.10))
  608 FORMAT (/,1X,T5,'hbar*omega for normal modes',//,1X,T5,'Bound',
     *   T20,1PE17.8, T40, 'Hartree', T50,1PE17.8, T70, 'eV', 
     *   /,1X,T20,1PE17.8, T40, 'kcal', T50,1PE17.8, T70,'cm**-1',
     *   /,1X,T5,'Unbound',T20,1PE17.8,T40,'i Hartree',T50,1PE17.8,
     *   T70,'i eV', /,1X,T20,1PE17.8, T40,'i kcal', T50,1PE17.8, 
     *   T70,'i cm**-1',//,1X,T5,'Morse parameters for bound motion')
  610 FORMAT (1X,T5,'DM', T33,1PE16.8, T50,'Hartree',T58,1PE16.8, 
     * T75,'eV', /,1X,T33,1PE16.8,T50,'kcal',
     * /,1X,T5,'Length parameter, AM', T33,1PE16.8, T50,'bohr**-1',
     * /,1X,T5,'KM=2.*SQRT(2.*MU*DM)/HBAR*AM',T33,1PE16.8,
     * /,1X,T5,'XM=1./KM',T33,1PE16.8)
  612 FORMAT (/,1X,T5,'WKB energy level for ground state',/,1X,T5,
     *   'E(N=  0)=',T15,1PE17.8, T33,'Hartree',T45,1PE17.8,T63,'eV', 
     * /,1X,T15,1PE17.8,T33,'kcal',T45,1PE17.8,T63,'cm**-1')
  614 FORMAT (/,1X,T5,'WKB energy level for excited state',/,1X,T5,
     * 'E(N=',I3, ')=',T15,1PE17.8,T33,'Hartree',T45,1PE17.8,T63,'eV', 
     * /,1X,T15,1PE17.8,T33,'kcal', T45,1PE17.8,T63, 'cm**-1')
  616 FORMAT (/,1X,T5,'Morse vibrational energy levels (relative to',
     *   ' bottom of vibrational well)'/,1X,T5,'(N = quantum number)',
     * //,1X,T4,'N',T11,'Hartree',T29,'eV',T47,'kcal',T65,'cm**-1')
  618 FORMAT (1X,T2,I5,T9,1PE16.8,T27,1PE16.8,T45,1PE16.8,T63,1PE16.8)
  620 FORMAT (1X,T5,'You have chosen IGS(6) = 1.  For a double ',
     * 'well potential',/,1X,T5,'the bend energy level is set to zero.')
  622 FORMAT (/, 1X, T5, 'Bending potential parameters at saddle point',
     *       //, 1X, T5, 'The potential for bend is approximated by',
     *        /, 1X, T10, 'V(PHI) = .5*FPHI*(PHI-PI)**2 + ',
     *                    '1./24.*APHI*(PHI-PI)**4',
     *        /, 1X, T12, 'FPHI = WB(BND)**2/GPHI',
     *       //, 1X, T69, 'Moment',
     *        /, 1X, T10, 'FPHI', T25, 'APHI', T36, 'hbar W(BND)', 
     *               T55, 'GPHI', T67, 'of inetria',                    TCA1097
     *        /, 1X, T10, '(kcal)', T24, '(kcal)', T38, '(cm**-1)', 
     *               T54, '(a.u.)', T69, '(a.u.)',
     *        /, 1X, T5, 1PE13.5, T20, 1PE13.5, T35, 1PE13.5,
     *               T50, 1PE13.5, T65, 1PE13.5)                        TCA1097
  624 FORMAT (/,1X,T5, 'Warning: the bending potential is a ',
     *                 'double well,',
     *        /,1X,T14, 'no harmonic approximation is tested',/)
  626 FORMAT (//,1X,T5, 'Coupled-mode bending energy levels at ',
     *                  'saddle point (kcal)',
     *        //,1X,T5, 'Harmonic', 
     *        (T24,1PE12.4,T37,1PE12.4,T50,1PE12.4,T62,1PE12.4))
627   FORMAT (1X,T5, 'Harmonic-quartic', 
     *        (T24,1PE12.4,T37,1PE12.4,T50,1PE12.4,T62,1PE12.4))
628   FORMAT (//, 1X,T5,'Single-mode bending energy levels at', 
     *   ' saddle point (kcal)',//,1X,T5, 'Harmonic',
     *   (T24,1PE12.4,T37,1PE12.4,T50,1PE12.4,T62,1PE12.4))
  630 FORMAT (/,1X,T5, 'Semiclassical solution,  range parameters ',
     *                 'not computed')
  632 FORMAT (/,1X,T5, 'Range parameters (radian**-1)',
     *        //,1X,T5, 'Harmonic',
     *   (T24,1PE12.4,T37,1PE12.4,T50,1PE12.4,T62,1PE12.4))
  634 FORMAT (/,1X,T5, 'Harmonic force constant < 0, cannot ',
     *                 'compute range parameters')
 6000 FORMAT (/,1X,T5,'Error: Root search to find the saddle point ',
     *                'has failed,',
     *        /,1X,T12,'the initial guess is probably bad.',
     *       //,1X,T5, 'R2', T70,'V(cal)'/)
 6001 FORMAT (/ 1X, F7.4, 2X, 11I10)
 6002 FORMAT (/ 5X, 'R1=', 2X, 11F10.4)
 6003 FORMAT (/,1X,T5,'Error: The roots of the force constant ',
     *                'matrix at the saddle point ',
     *        /,1X,T12,'are not opposite in sign, or one is zero.')
      END
