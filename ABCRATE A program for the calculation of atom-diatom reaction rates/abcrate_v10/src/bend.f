!!!*************************************************************
! 文件/File: bend.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: bend.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:40
!*************************************************************

C
      SUBROUTINE BEND (X, Y, FB, AB, GB, XMOM)
C
C     BEND - compute bending potential parameters
C
C  Called by:
C     DATAIN - read in data and set up constants
C     PARAM  - compute parameters for bound motion along reaction path
C
C  Calls:
C     BENPOT - evaluates bending potential
C     BENMIN - find minimum of bending potential
C     LIN2   - solve two simultaneous linear equations
C     PEF    - potential subroutine
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LCONV
      DIMENSION VBEND(2)
      COMMON /BENDOP/ THETA(2), A11, A12, A21, A22, IBOPT
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
C
C  Compute potential on MEP
C     WRITE (6, 6600) X,Y,THETA, A11, A12, A21, A22, IBOPT
C6600 FORMAT (' ENTER BEND: X,Y,THETA,A11,A12,A21,A22,IBOPT='/
C    *   1P8E15.7, I5)
      R2 = Y/CM2
      R1 = X - CM1*R2
      R3 = R1+R2
      CALL PEF (R1, R2, R3, V, D1, D2, D3, 1)                           GCL0695
      R1S = R1*R1
      R2S = R2*R2
      TR12 = 2.D0*R1*R2
      VCOL=V
      IF (IBOPT .GE. 3) THEN
C  fit 3; fit bend parameters to bending potential evaluated at two
C     angles THETA(1-2)
         DO 10 I=1,2
            PHI = THETA(I)
            CALL BENPOT (PHI, R1S, R2S, TR12)
            VBEND(I)=V-VCOL
   10    CONTINUE
         CALL LIN2 (A11, A12, A21, A22, VBEND(1), VBEND(2), FB, AB)
      ELSE
C  fits 1 and 2; fit harmonic potential using first derivative
         T1 = R1*R2/R3
         R3S = R3
         DD1 = D3
         FB = -T1*DD1
         IF (FB .LT. 0.D0) THEN
C  FB is less than zero, fit to double well; first search for minimum in
C     bending potential
            AB = 0.D0                                                   GCL1092
            CALL BENMIN (VMIN, PHIMN, IERR, R1S, R2S, TR12)
            IF (IERR .EQ. 0) THEN
               VB = VMIN - VCOL
               IF (IBOPT .EQ. 1) THEN
C  Fit 2a, fit bend to give VB and correct tp for E=VCOL; first find tp
C     for E=VCOL
                  SGN = -1.0D0
                  PHI = PHIMN
                  DELPHI = MIN(PI/18.D0, 0.5D0*PHI)
                  AB = 0.D0                                             GCL1092
                  LCONV = .FALSE.
C                 WRITE (6, 6620)
C6620 FORMAT (' FIT 2A, FIND V(PHI)=VCOL'/ 8X, 'PHI', 12X, 'V', 14X,
C    *   'VCOL')
   20             CONTINUE
C  PHI is restricted to be between 0 and pi, therefore want tp for
C     angle lower than pi.
                     PHI = PHI - DELPHI
                     IF (PHI .LT. 0.0D0) THEN
                        WRITE (6, 601)
                        GO TO 30
                     END IF
                     CALL BENPOT (PHI, R1S, R2S, TR12)
C                    WRITE (6, 6621) PHI, V, VCOL
C6621 FORMAT (1X, 1P3E17.9)
                     IF (SGN*(V-VCOL) .GT. 0.D0) GO TO 20
                     T = V - VCOL
                     IF (ABS(VCOL) .GT. 1.D-8) T = T/VCOL
                     LCONV = ABS(T) .LE. 1.D-8 .OR.
     *                  ABS(DELPHI) .LT. 1.D-10
                     SGN = - SGN
                     DELPHI = -0.5D0*DELPHI
                  IF (.NOT. LCONV) GO TO 20
   30             CONTINUE
C  If couldn't find turning point FB is left negative and AB is left
C     zero.
                  IF (LCONV) THEN
C  Turning point found, fit parameters
                     T = PHI - PI
                     IF (ABS(T) .GT. 1.D-8) THEN
                        T2 = T*T
                        T4 = T2*T2
                        FB = 8.D0*VB/T2
                        AB = -96.D0*VB/T4
                     END IF
                  END IF
               END IF
C  If AB is zero either IBOPT=2 or fit 2a failed.
               IF (AB .EQ. 0.D0) THEN
C  Fit 2b,  fit bend to give potential and location at minimum
C                WRITE (6, 6630) VB, PHIMN
C6630 FORMAT (' FIT 2B, VB, PHIMN=', 1P2E17.9)
                  T = PHIMN - PI
                  IF (ABS(T) .GT. 1.D-8) THEN
                     T2 = T*T
                     T4 = T2*T2
                     FB = 4.0D0*VB/T2
                     AB = -24.D0*VB/T4
                  END IF
               END IF
            END IF
         ELSE
C  Fit 1, FB > 0, fit AB to give correct 4th derivative w.r.t. angle at
C     minimum.  4th derivative wrt angle is proportional to the 2nd
C     derivative wrt R3; found by numerical first derivative of analytic
C     first derivative.
            DEL = 1.D-1                                                 TCA1195
            DD2 = 0.D0                                                  GCL1092
            N = 0
            CALL PEF (R1,R2,R3,V,D1,D2,D3,1)                            TCA1195
            DX0=D3                                                      TCA1295
C  Loop over decreasing step size until converged
   40       CONTINUE
               N = N + 1
               D2OLD = DD2
               DEL = DEL*.1D0
               R3 = R3S - DEL
               CALL PEF (R1, R2, R3, V, D1, D2, D3, 1)                  GCL0893
               DX1=D3                                                   TCA1295
               R3 = R3S - 2.D0*DEL                                      TCA1295
               CALL PEF(R1,R2,R3,V,D1,D2,D3,1)                          TCA1295
               DX2=D3                                                   TCA1295
               DD2 = 0.5D0*(DX2-4.D0*DX1+3.D0*DX0)/DEL                  TCA1295
               T = DD2-D2OLD
               IF (ABS(DD2) .GT. 1.D-8) T = T/DD2
               LCONV = ABS(T) .LT. 1.D-6
            IF (.NOT.LCONV .AND. N .LT. 7) GO TO 40                     TCA1195
            IF (.NOT. LCONV) WRITE (6,600) DD2, D2OLD
            T2 = 3.D0*T1
            AB = (T2*DD2 + (1.D0-T2)*DD1)*T1
         END IF
      END IF
C  Moment of inertia
      R1I = 1.D0/R1
      R2I = 1.D0/R2
      GB = R1I*R1I/XMA + R2I*R2I/XMC + ((R1I+R2I)**2)/XMB
      XMOM = XMA*XMB*XMC*R1*R2*R1*R2*GB/(XMA+XMB+XMC)
C      WRITE (6, 6650) FB, AB
C6650  FORMAT (' EXIT BEND: FB, AB=', 1P2E15.7)
      RETURN
C
600   FORMAT(/,2X,T5,'Warning: Quartic bend not converged, the ',
     *               'last two iterations are: ',2(1PE13.5,1X)) 
601   FORMAT(/,2X,T5,'Warning: In fit of the bend pot., could ',
     *               'not find turning pt. to be used')
      END

