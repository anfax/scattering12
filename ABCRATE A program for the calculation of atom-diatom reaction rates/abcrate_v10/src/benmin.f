!!!*************************************************************
! 文件/File: benmin.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: benmin.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE BENMIN (VMIN, PHIMN, IERR, R1S, R2S, TR12)
C
C     BENMIN - find minimum of bending potential
C
C  Called by:
C     BEND - compute bending potential parameters
C
C  Calls:
C     BENPOT - evaluates bending potential
C     QUADFT - quadratic fit of three points
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LCONV
      DIMENSION PHIQ(3), VQ(3), B(3)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
C
      VMIN = V
C     WRITE (6, 6601) PI, V
C6601 FORMAT (' FIND MINIMUM IN V(PHI)'/ 1X, 'IC', 4X, 'PHI', 14X, 'V',
C    *   16X, 'B', 44X, 'PHIL', 13X, 'PHIR'/ 3X, 1P2E17.9)
      DELPHI = PI/18.D0                                                 GCL1092
      PHI = PI
      PHIQ(1) = PI
      PHIQ(2) = PI
      VQ(1) = V
      VQ(2) = V
      PHIMN = PI
      IC = 1
      IERR = 1
C  First scan from collinear (DEL=PI) to lower angles to see if the
C     potential drops as the angle decreases.
   10 CONTINUE
         IC = IC + 1
         PHI = PHI - DELPHI
C  If the potential was flat or continually drop from the collinear
C     geometry then assume a minimum doesn't exits
         IF (PHI .LE. 0.D0) RETURN
         CALL BENPOT (PHI, R1S, R2S, TR12)
C        WRITE (6, 6602) PHI, V
C6602 FORMAT (3X, 1P2E17.9)
         IF (V .GT. VMIN) GO TO 20
         PHIQ(1) = PHIQ(2)
         PHIQ(2) = PHI
         PHIMN = PHI
         VQ(1) = VQ(2)
         VQ(2) = V
         VMIN = V
         GO TO 10
   20 CONTINUE
C  PHIQ(1) and PHIQ(2) are bounds on location of the minimum
      IC = MIN(3, IC)
      PHIQ(IC) = PHI
      VQ(IC) = V
C  need three points in the root search below, if only have two get the
C     third
      IF (IC .LT. 3) THEN
         PHI = 0.5D0 * (PHIQ(1) + PHIQ(2))
         CALL BENPOT (PHI, R1S, R2S, TR12)
C        WRITE (6, 6602) PHI, V
         PHIQ(3) = PHI
         VQ(3) = V
         PHIL = PHIQ(2)
         PHIR = PHIQ(1)
         IF (V .GT. VMIN) THEN
            PHIL = PHIQ(3)
         ELSE
            PHIMN = PHIQ(3)
            VMIN = V
         END IF
      ELSE
         PHIL = PHIQ(3)
         PHIR = PHIQ(1)
      END IF
C  start root search for phi
      IC = 0
      LCONV = .FALSE.
   30 CONTINUE
         IC = IC + 1
         PHIOLD = PHI
         PHIMNO = PHIMN
         VMINO = VMIN
C   try quadratic fit
         CALL QUADFT (PHIQ, VQ, B)
C        WRITE (6, 6604) IC, B, PHIL, PHIR
C6604 FORMAT (1X, I2, 34X, 1P3E15.7, 1P2E17.9)
         VV = 0.0D0
         IF (B(3) .GT. 0.D0) THEN
C   extrema is a min, use quadratic fit unless out of range
            PHI = -0.5D0 * B(2) / B(3)
            IF (PHI .LE. PHIL) THEN
               PHI = 0.5D0 * (PHIL + PHIMN)
            ELSE IF (PHI .GE. PHIR) THEN
               PHI = 0.5D0 * (PHIMN + PHIR)
            ELSE
               VV = B(1) + PHI * (B(2) + B(3) * PHI)
            END IF
         ELSE
C   otherwise average left and right bounds
            PHI = 0.5D0 * (PHIL + PHIR)
         END IF
C
         CALL BENPOT (PHI, R1S, R2S, TR12)
C        WRITE (6, 6602) PHI, V
         IF (V .GT. VMIN) THEN
C   V > VMIN, reset left or right bound
            IF (PHI .LT. PHIMN) THEN
               PHIL = PHI
            ELSE
               PHIR = PHI
            END IF
         ELSE
            VMIN = V
            PHIMN = PHI
            T = VMINO - V
            IF (ABS(VMINO) .GT. 1.D-8) T = T/VMINO
            LCONV = ABS(T) .LE. 1.D-8
C   V< VMIN, reset left or right bound
            IF(PHIMNO .LT. PHIMN) THEN
               PHIL = PHIMNO
            ELSE
               PHIR = PHIMNO
            END IF
         END IF
C   update PHIQ and VQ with most current guesses
         I = 3
         DO 40 II = 1,2
            I = I - 1
            PHIQ(I+1) = PHIQ(I)
            VQ(I+1) = VQ(I)
   40    CONTINUE
         PHIQ(1) = PHI
         VQ(1) = V
C
C   check for convergence in PHI or V
         IF (.NOT.LCONV) THEN
            LCONV = ABS(PHI - PHIOLD) .LT. 1.D-6
            IF (.NOT.LCONV) THEN
               T = V - VV
               IF (ABS(V) .GT. 1.D-8) T = T/V
               LCONV = ABS(T) .LE. 1.D-8
            END IF
         END IF
      IF (.NOT. LCONV .AND. IC .LT. 30) GO TO 30
C     WRITE (6, 6610) PHI, PHIMN, V, VMIN
C6610 FORMAT (' PHI, PHIMN, V, VMIN=', 1P4E15.7)
      IF (LCONV) IERR = 0
      RETURN
      END
