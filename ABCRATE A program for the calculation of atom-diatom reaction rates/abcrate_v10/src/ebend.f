!!!*************************************************************
! 文件/File: ebend.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: ebend.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      FUNCTION EBEND (N, FB, AB, GB)
C
C     EBEND  - compute bending energy levels
C
C  Called by:
C     ADIAB  - compute adiabatic potential curves
C     GTST   - compute free energies, CVT and ICVT rates
C     PRTQRT - compute partition function of harm-quartic potential
C     RPHSUM - summarize reaction path info
C     SADDLE - find saddle point(s) and do normal mode analysis
C     TABL21 - print out table of GTS info
C     THRCOR - compute threshold corrections for ICVT
C     VBEND  - compute bending energy levels
C     VSPLN2 - spline fit of adiabatic potential in product channel
C
C  Calls:
C     DBWL   - compute semiclassical eigenvalue of harm.-quart. with
C              double well
C     HQSC   - compute semiclassical eigenvalues of harm.-quart.
C              potential without double well.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LFLG
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      LOGICAL LMAX
      COMMON /EBND1/   LMAX
      COMMON /EBND2/  ALF, ALFBAR
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /MORBC/  DBC, XKBC, AMBC
      SAVE C3                                                           TCA1097
      DATA C3 /0.333333333333D0/
      SAVE ITYP,XMU,V0,WB,CC,AB24,RTC3,EMAX,ELAST,ALFSQ                 TCA1097
C
      EBEND = 0.D0                                                      GCL1092
      ALF = 0.D0                                                        GCL1092
      ALFBAR = 0.D0                                                     GCL1092
      ALFSQ = 0.0D0
      XLAMSQ = 1.0D0
      LMAX = .FALSE.
      LFLG = .TRUE.
      IF (LGS(4)) THEN
C  Harmonic-quartic
         E = 0.D0                                                       GCL1092
         ELAST = -1.D35
         V0 = 0.D0                                                      GCL1092
         IF (FB .LE. 0.D0 .AND. AB .LE. 0.D0) THEN                      GCL1092
            LMAX = .TRUE.
            RETURN
         END IF
         IF (LGS(8) .OR. FB .LE. 0.D0) THEN
C  set up for semiclassical solution
            XMU = 1.0D0/GB
            A = 0.5D0*ABS(FB)
            B = AB/24.0D0
            EMAX = DBC
            EM = 1.D-5
            IF (AB .LT. 0.D0) EMAX = -0.25D0*A*A/B
            IF (FB .LT. 0.D0) V0 = .25D0*A*A/B
            IF (AB .LE. 0.D0) ITYP=1
            IF (AB .GT. 0.D0 .AND. FB .GE. 0.D0) ITYP=2                 GCL1092
            IF (AB .GT. 0.D0 .AND. FB .LT. 0.D0) ITYP=3                 GCL1092
         ELSE
C  set up for perturbation-variation solution for FB>0
            RTC3 = SQRT(C3)
            WB = SQRT(FB*GB)
            ALFSQ = SQRT(FB/GB)
            E0 = 0.D0                                                   GCL1092
            EMAX=DBC
            CC = GB/WB
            CC = 0.75D0*CC*CC
            AB24 = AB/24.D0                                             GCL1092
            EMAX = DBC
            T = -6.D0*FB/AB
            IF (T .GE. 0.D0) EMAX = 0.5D0*T*(FB+2.D0*AB24*T)
         END IF
      ELSE
C  harmonic
         IF (FB .LE. 0.D0) THEN
            LMAX = .TRUE.
            RETURN
         ELSE
            WB = SQRT(FB*GB)
         END IF
      END IF
      GO TO 1
C
      ENTRY EBN(N, FB, AB, GB)
C
      LFLG = .FALSE.
    1 CONTINUE
      IF (LGS(4)) THEN
C  harmonic-quartic potential
         IF (LGS(8) .OR. FB .LE .0.D0) THEN
            IF (ITYP .EQ. 1) THEN
C   semiclassical solution for AB<0
               CALL HQSC(XMU, A, B, N, E)                               GCL1096
            ELSE
C   semiclassical solution for AB>0
               IF (LFLG .OR. (-1)**N .GT. 0.D0) THEN
                  NQ = N/2 + 1
                  CALL DBWL (XMU, A, B, V0, NQ, EM, EP, ITYP)
               END IF
               IF ((-1)**N .GT. 0.D0) THEN
                  E = EM - V0
               ELSE
                  E = EP - V0
               END IF
            END IF
         ELSE
C   perturbation-variation soltuion for FB>0
            FN = DBLE(N)                                                GCL1092
            E0 = (FN+0.5D0)*WB
            E = E0
            XLAMSQ = 0.D0                                               GCL1092
            E1 = CC*(2.D0*FN*(FN+1.D0) + 1.D0)
            C1 = 2.D0*AB24*E1/E0
            IF (C1 .NE. 0.D0) THEN
               C2 = C1*C1 - 1.D0/27.0D0
               IF (C2 .LT. 0.0D0) THEN
                  PHI = C3*ACOS(3.D0*C1/RTC3)
                  XLAMSQ = 2.0D0*RTC3*COS(PHI)
               ELSE
                  T2 = SQRT(C2)
                  T = C1 + T2
                  T1 = C1 - T2
                  XLAMSQ = SIGN(ABS(T)**C3, T) + SIGN(ABS(T1)**C3, T1)
               END IF
               T = 1.D0/XLAMSQ
               E = 0.5D0*(XLAMSQ+T)*E0 + AB24*E1*T*T
            END IF
         END IF
         EBEND = E
         EBN = E
         IF (E .GT. EMAX .OR. E .LT. ELAST .OR. E .LT. -V0) THEN
            WRITE(6, 600) N, E, EMAX, ELAST, -V0
            IF (E .LT. -V0) EBEND = -V0
            LMAX = .TRUE.
         ELSE
            ELAST = E
            IF (FB .GT. 0.D0) THEN
               ALF = SQRT(ALFSQ)
               ALFBAR = SQRT(XLAMSQ*ALFSQ)
            END IF
         END IF
      ELSE
C  harmonic
         E = (DBLE(N) + 0 .5D0)*WB                                      GCL1092
         EBN = E
         EBEND = E
      END IF
      RETURN
600   FORMAT(/,2X,T5,'Warning: In EBEND N = ',I4,' E = ',1PE13.5,
     *               ' EMAX = ',1PE13.5,' previous E = ',1PE13.5,
     *               ' V0 = ',1PE13.5)
      END
