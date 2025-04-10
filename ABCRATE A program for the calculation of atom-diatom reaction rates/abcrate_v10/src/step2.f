!!!*************************************************************
! 文件/File: step2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: step2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE STEP2 (X, Y, DX, DY, STEP, V0)
C
C     STEP2  - follow gradient by choosing step to minimize potential
C              along an arc
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Calls:
C     DERIV  - compute derivatives of potential w.r.t. x,y coordinates
C     QUADFT - quadratic fit of three points
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LCONV
      DIMENSION PHIQ(3), VQ(3), B(3)
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
C
      XZ = X
      YZ = Y
      PHI0 = ATAN2 (DY, DX)
      PHI = PHI0
      CALL ST2DRV (DX, DY, XZ, YZ, XD, YD, STEP, PHI)
      PHIMIN = PHI
      XDMIN = XD
      YDMIN = YD
      VMIN = V
      DELPHI = 0.25D0
      PHIQ(1) = PHI
      PHIQ(2) = PHI
      VQ(1) = V
      VQ(2) = V
      IC = 1
C     WRITE (91, 9100) PHI0, VMIN
C9100 FORMAT (' STEP2: PHI0,VMIN=',1P2E13.5)
   10 CONTINUE
         IC = IC + 1
         PHI = PHI + DELPHI
         CALL ST2DRV (DX, DY, XZ, YZ, XD, YD, STEP, PHI)
C        WRITE (91, 9101) PHI, V
C9101    FORMAT (8X, 'PHI ,V   =', 1P2E13.5)
         IF (V .GT. VMIN) GO TO 20
         PHIQ(1) = PHIQ(2)
         PHIQ(2) = PHI
         PHIMIN = PHI
         XDMIN = XD
         YDMIN = YD
         VMIN = V
         VQ(1) = VQ(2)
         VQ(2) = V
         VMIN = V
      GO TO 10
   20 CONTINUE
      IC = MIN(3, IC)
      PHIQ(IC) = PHI
      VQ(IC) = V
C     WRITE (91, 9102) IC, (PHIQ(I), VQ(I), I=1,IC)
C9102 FORMAT (1X, I5, ' POINTS'/ (8X, 'PHI ,V   =',1P2E13.5))
      IF (IC .LT. 3) THEN
         I = 3
         DO 30 II = 1,2
            I = I - 1
            PHIQ(I+1) = PHIQ(I)
            VQ(I+1) = VQ(I)
   30    CONTINUE
         PHI = PHI0
         DELPHI = -DELPHI
   40    CONTINUE
            PHI = PHI + DELPHI
            CALL ST2DRV (DX, DY, XZ, YZ, XD, YD, STEP, PHI)
C           WRITE (91, 9101) PHI, V
            IF (V .GT. VMIN) GO TO 50
            PHIQ(3) = PHIQ(2)
            PHIQ(2) = PHI
            PHIMIN = PHI
            XDMIN = XD
            YDMIN = YD
            VMIN = V
            VQ(3) = VQ(2)
            VQ(2) = V
            VMIN = V
         GO TO 40
   50    CONTINUE
         PHIQ(1) = PHI
         VQ(1) = V
      END IF
      PHIL = PHIQ(1)
      PHIR = PHIQ(3)
C  start root search for phi
      IC = 0
      LCONV = .FALSE.
   60 CONTINUE
         IC = IC + 1
         PHIOLD = PHI
         PHIMNO = PHIMIN
         VMINO = VMIN
C   try quadratic fit
         CALL QUADFT (PHIQ, VQ, B)
         IF (B(3) .GT. 0.D0) THEN
C   extrema is a min, use quadratic fit unless out of range
            PHI = -0.5D0 * B(2) / B(3)
            IF (PHI .LE. PHIL .OR. PHI .GE. PHIR) THEN
               PHI = 0.5D0 * (PHIL + PHIR)
            ELSE
               VV = B(1) + PHI * (B(2) + B(3) * PHI)
            END IF
         ELSE
            PHI = 0.5D0 * (PHIL + PHIR)
         END IF
C
         CALL ST2DRV (DX, DY, XZ, YZ, XD, YD, STEP, PHI)
C        WRITE (91, 9104) IC, PHI, PHIL, PHIR, V, B
C9104    FORMAT (' IC,PHI,PHIL,PHIR,V,B=', I5, 1P7E13.5)
         IF (V .GT. VMIN) THEN
            IF (PHI .LT. PHIMIN) THEN
               PHIL = PHI
            ELSE
               PHIR = PHI
            END IF
         ELSE
            IF(PHIMNO .LT. PHIMIN) THEN
               PHIL = PHIMNO
            ELSE
               PHIR = PHIMNO
            END IF
            PHIMIN = PHI
            VMIN = V
            XDMIN = XD
            YDMIN = YD
         END IF
         I = 3
         DO 70 II = 1,2
            I = I - 1
            PHIQ(I+1) = PHIQ(I)
            VQ(I+1) = VQ(I)
   70    CONTINUE
         PHIQ(1) = PHI
         VQ(1) = V
C
C   check for convergence in PHI or V
         LCONV = ABS(PHI - PHIOLD) .LT. 1.D-6
         IF (.NOT.LCONV) THEN
            LCONV = ABS(1.D0 - VV / V) .LE. 1.D-8 .OR.
     *         ABS(1.D0 - V / VMINO) .LE. 1.D-8
         END IF
      IF (.NOT.LCONV .AND. IC .LT. 30) GO TO 60
   80 CONTINUE
      CALL DERIV (XDMIN, YDMIN, DX, DY)
      DMAG = SQRT(DX*DX+DY*DY)
      DX = DX/DMAG
      DY = DY/DMAG
      X = XDMIN
      Y = YDMIN
      V0 = VMIN
      RETURN
      END
