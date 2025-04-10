!!!*************************************************************
! 文件/File: geom.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: geom.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE GEOM (LP, AA, AB, IND, VINF, RASY, V0, DM, OM, XKM, AM)
C
C     GEOM   - find equilibrium geometry for reactants and products
C
C  Called by:
C     DATAIN - read in data and set up constants
C
C  Calls:
C     PEF    - evaluate potential
C     D2VDU2 - compute second derivative along u coordinate
C     VMIN   - find minimum energy along u coordinate
C     WKBSET - set up grid of WKB energy levels
C     WKBINT - extrapolate WKB energy levels from grid
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*2 AA, AB                                                GCL1096
      LOGICAL LP
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      LOGICAL LGS2(10)
      COMMON /LOGIC2/ LGS2
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
      COMMON /STATE/  TNP1, LSTATE, NSTATE
C
      IF (IND .LE. 0) THEN
C  reactants
         IS = -1
         S = -1000.D0                                                   GCL1092
         IF (LP) WRITE (6, 600)
         DX = 1.D0                                                      GCL1092
         DY = 0.D0                                                      GCL1092
         R1 = 1000.D0                                                   GCL1092
         R2 = RASY
         RSAVE = R2
      ELSE
C  products
         IS = -2
         S = 1000.D0                                                    GCL1092
         T = SQRT(CM1*CM1 + CM2*CM2)
         DX = CM1/T
         DY = CM2/T
         IF (LP) WRITE (6, 602)
         R2 = 1000.D0                                                   GCL1092
         R1 = RASY
         RSAVE = R1
      END IF
      IC = 0
C  Loop over asymptotic distance
   10 CONTINUE
         IF (IC .GT. 10) THEN
            WRITE (6, 6000)
            R = RSAVE - 1.0D0
            IF (R .LT. 0.D0) R= 0.D0                                    GCL1092
            DO 20 I = 1,20
               R = R + .1D0
               IF (IND .GT. 0) R1 = R
               IF (IND .LE. 0) R2 = R
               R3 = R1 + R2
               CALL PEF (R1, R2, R3, V, D1, D2, D3, 0)                  GCL0893
               WRITE (6, 6001) R1, R2, V
   20       CONTINUE
            GO TO 30
         END IF
         IF (IND .GT. 0) R2 = R2 + 100.D0                               GCL1092
         IF (IND .LE. 0) R1 = R1 + 100.D0                               GCL1092
         Y = CM2*R2
         X = R1 + CM1*R2
C  find minimum
         CALL VMIN (X, Y, DX, DY, 0.1D0, V0, .FALSE., IERR)
         IF (IERR .NE. 0) THEN
            WRITE (6, 6002)
            STOP 'GEOM 1'
         END IF
         IF (IND .GT. 0) RASY = X-CM1*Y/CM2
         IF (IND .LE. 0) RASY = Y/CM2
         IC = IC + 1
         IF (IC .GT. 1) THEN
C  check for convergence
            T1 = RASY-ROLD
            IF (RASY .NE. 0.D0) T1 = T1/RASY
            IF (ABS(T1) .LT. 1.D-8) THEN
               T1 = VOLD-V0
               IF (ABS(V0) .GT. 1.D-6) T1=T1/V0
               IF (ABS(T1) .LT. 1.D-10) GO TO 30
            END IF
         END IF
         ROLD = RASY
         VOLD = V0
      GO TO 10
   30 CONTINUE
C  get Morse parameters
      UX = -DY
      UY = DX
      CALL D2VDU2 (X, Y, UX, UY, DD2)
      DM = VINF - V0
      OM = SQRT(DD2/XMU)
      XKM = 4.D0*DM/OM
      XE = 0.D0                                                         GCL0795
      IF (IND .GT. 0) XMX = XMP
      IF (IND .LE. 0) XMX = XM
      AM = 2.D0*SQRT(2.D0*XMX*DM)/XKM
      IF (LGS(7) .AND. .NOT.LGS2(5)) THEN
C  compute WKB energy level
         E = 0.5D0*OM*(1.D0-0.5D0*XE)
         CALL WKBSET (IS, S, 0, E, UL, UG, PER)
         IF (LSTATE .NE. 0 .AND. NSTATE .NE. 0) THEN
            E = 0.5D0*TNP1*OM*(1.D0-0.5D0*TNP1*XE)
C  if state-selected and NSTATE>0 compute ground-state WKB energy
            CALL WKBSET (IS, S, NSTATE, E, UL, UG, PER)
         END IF
      END IF
      IF (.NOT.LP) RETURN
C  print out
      T1 = OM*EAU
      T2 = OM*CKCAL
      T3 = OM*CCM
      WRITE (6, 604) AA, AB, RASY, AA, AB, OM, T1, T2, T3, AA, AB
      T1 = DM*EAU
      T2 = DM*CKCAL
      XE = 1.D0/XKM
      WRITE (6, 606) DM, T1, T2, AM, XKM, XE
C  diatomic vibrational energy levels
      IF (LGS(7)) THEN
         CALL WKBINT (IS, S, 0, E, UL, UG, PER)
         WRITE (6, 608) E, E*EAU, E*CKCAL, E*CCM
         IF (LSTATE .NE. 0 .AND. NSTATE .NE. 0) THEN
            CALL WKBINT (IS, S, NSTATE, E, UL, UG, PER)
            WRITE (6, 610) NSTATE, E, E*EAU, E*CKCAL, E*CCM
         END IF
      END IF
      WRITE (6, 612)
      NMAX = INT(.5D0*(XKM-1)) + 1
      TIXK = 2.D0/XKM
      T = -.5D0*TIXK
      DO 40 N = 1, NMAX
         NQM = N-1
         T = T + TIXK
         E = DM*T*(2.D0-T)
         T1 = E*EAU
         T2 = E*CKCAL
         T3 = E*CCM
         WRITE (6, 614) NQM, E, T1, T2, T3
   40 CONTINUE
      RETURN
  600 FORMAT (/,1X,28('*'),2X,'Reactant properties',1X, 28('*') )
  602 FORMAT (/,1X,29('*'),1X,'Product properites',1X, 29('*') )
  604 FORMAT (/,1X,T5,2A2, ' equilibrium distance',T33,F12.8,T50,'bohr',
     *        /,1X,T5,'hbar*omega of ', 2A2, ' vibration',T33,1PE16.8, 
     * T50,'Hartree', T58,1PE16.8, T75,'eV', 
     * /,1X,T33, 1PE16.8, T50, 'kcal', T58,1PE16.8, T75,'cm**-1',
     * //,1X,T5,'Morse parameters, ', 2A2, ' vibration')
  606 FORMAT (1X,T5,'DM', T33,1PE16.8, T50,'Hartree',T58,1PE16.8, 
     * T75,'eV', /,1X,T33,1PE16.8,T50,'kcal',
     * /,1X,T5,'Length parameter, AM', T33,1PE16.8, T50,'bohr**-1',
     * /,1X,T5,'KM=2.*SQRT(2.*MU*DM)/HBAR*AM',T33,1PE16.8,
     * /,1X,T5,'XM=1./KM',T33,1PE16.8)
  608 FORMAT (/,1X,T5,'WKB energy level for ground state',/,1X,T5,
     *   'E(N=  0)=',T15,1PE17.8, T33,'Hartree',T45,1PE17.8,T63,'eV', 
     * /,1X,T15,1PE17.8,T33,'kcal',T45,1PE17.8,T63,'cm**-1')
  610 FORMAT (/,1X,T5,'WKB energy level for excited state',/,1X,T5,
     * 'E(N=',I3, ')=',T15,1PE17.8,T33,'Hartree',T45,1PE17.8,T63,'eV', 
     * /,1X,T15,1PE17.8,T33,'kcal', T45,1PE17.8,T63, 'cm**-1')
  612 FORMAT (/,1X,T5,'Morse vibrational energy levels (relative to',
     *   ' bottom of vibrational well)'/,1X,T5,'(N = quantum number)',
     * //,1X,T4,'N',T14,'Hartree',T34,'eV',T51,'kcal',T68,'cm**-1')
  614 FORMAT (1X,T2,I5,T9,1PE16.8,T27,1PE16.8,T45,1PE16.8,T63,1PE16.8)
 6000 FORMAT (/,1X,T5,'The calculation of asymptotic geometry does ',
     * 'not converge with respect to',/,1X,T5,'moving further out ',
     * 'into the asymptotic region',//,1X,T10,'R1',T25,'R2',T35,'V')
 6001 FORMAT (1X,T5,F12.6,T20,F12.6,T35,1PE14.6)
6002  FORMAT(/,2X,T5,'Error: In GEOM there is a problem with the call ',
     *               'to the VMIN subprogram.')
      END
