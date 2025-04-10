!!!*************************************************************
! 文件/File: hqsc.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: hqsc.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE HQSC (XMU, A, B, NQ, E)
C
C     HQSC   - compute semiclassical eigenvalues of harm.-quart.
C              potential.  Assumes A>0, B<0.
C
C  Called by:
C     EBEND  - compute bending energy levels
C
C  Calls:
C     ELLIP  - elliptic integrals of first and second kind
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      SAVE EPS                                                          TCA1097
      DATA EPS/1.D-8/
C
      XN = DBLE(NQ) + 0.5D0                                             GCL1092
      IF (B .EQ. 0.D0) THEN
C  harmonic for B=0
         E = XN*SQRT(2.0D0*A/XMU)
      ELSE
C  harmonic-quartic
         BP = -B
         T1 = 0.5D0*A
         V0 = T1*T1/BP
         T1 = SQRT(XMU*A)
         THMX = 2.0D0*T1*A/(3.0D0*PI*BP)
         DTHC = 2.0D0*T1/(PI*A)
         FN1 = THMX - XN
         IF (FN1 .LT. 0.0D0) THEN
C  Quantum number is above the maximum; set energy to V0 and return
            E = 1.00001D0*V0
         ELSE
C  Linear approximation to energy root
            E2 = XN*V0/THMX
            EL = 0.0D0
            EG = V0
C  Root search for energy level
   10       CONTINUE
C     Evaluate phase integral
               T1 = SQRT ((V0-E)/V0)
               T2 = 1.0D0+T1
               RTT2 = SQRT(T2)
               RS = (1.0D0-T1)/T2
               CALL ELLIP (RS, EL1, EL2)
               TH = THMX*RTT2*(EL2 - T1*EL1)
               FN2 = TH - XN
               IF (ABS(FN2) .GT. EPS) THEN
C     Not converged, get next approximation to the energy
                   IF (FN2 .GT. 0.0D0) EG = MIN(E2,EG)
                   IF (FN2 .LT. 0.0D0) EL = MAX(E2,EL)
                   E = 0.0D0
                   IF (ABS(FN2) .LT. ABS(FN1)) THEN
C     Take Newton-Raphson step; need derivative of phase integral
                      DTH = DTHC*EL1/RTT2
                      IF (DTH .NE. 0.0D0) E = E2 - FN2/DTH
                   END IF
C     If we can't take a Newton-Raphson step take average of the minimum
C        and maximum energies.
                   IF (E .LE. EL .OR. E .GE. EG) E = 0.5D0*(EL + EG)
                   E2 = E
                   FN1 = FN2
C  Loop back over interations in the root search
                   GO TO 10
                END IF
         END IF
      END IF
      RETURN
      END
