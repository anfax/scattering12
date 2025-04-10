!!!*************************************************************
! 文件/File: vsplin.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vsplin.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE VSPLIN (NSMAX, SS, SMAX, VMAX, VAD)
C
C     VSPLIN - spline fit of adiabatic potential
C
C  Called by:
C     KAPVA  - compute kappas
C
C  Calls:
C     SPL1D1 - spline fit
C     SPL1B1 - spline fit
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      DIMENSION IOP(2), SS(NSDM10), VAD(NSDM)                           TCA0197
      COMMON /SPLNVA/ SV(NSDM1), VV(NSDM1), AV(NSDM1), BV(NSDM1),
     *CV(NSDM1), DV(NSDM1), SCR(NSDM1), NS
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
C
      NS = 0
      NSLAST = NSMAX
      IF (LSYM) NSLAST = NMID
      IORG = 1
C  IOP(1)=5 specifies that a 4-pt. Lagrance interpolation will be used
C     to determine derivative at left edge of grid.
      IOP(1) = 5
      SCR(1) = 0.0D0
C  load SV and VV arrays for spline fitting with adiabatic potential
      IS = 0
      DO 10 I = 1,NSLAST
         IS = IS + 1
         SV(IS) = SS(IS)
         VV(IS) = VAD(IS)
         IF (SS(IS) .GE. SMAX) GO TO 20
   10 CONTINUE
   20 CONTINUE
C  check if maximum should be added to the arrays
      IF (IS .LE. 1 .OR. (.NOT.LSYM .AND. IS .GE. NSLAST-1)) THEN       TCA1296
C  maximum at an edge, skip adding it
         IS0 = 0
         IORG = 1
      ELSE
         IS0 = 1
C  ISS is first point past (or equal to) the location of the maximum
         ISS = IS
         IF (ABS(SS(IS-1)-SMAX) .LT.  1.D-4) THEN
C  maximum is within 1.e-4 of the ISS-1 grid point, therefore don't add
C     an extra point, just replace the point at ISS-1 with the values at
C     the maximum; IS is now set to the previous grid point for spline
C     fit
            IS0 = 0
            ISS = IS - 1                                                TCA1296
         ELSE IF (ABS(SS(IS)-SMAX) .LT.  1.D-4) THEN
C  maximum is within 1.e-4 of the grid point, therefore don't add an
C     extra point, just replace the point at IS; ISS is now set to the
C     next grid point
            IS0 = 0
            ISS = IS + 1
         END IF
C  replace IS location with maximum
         SV(IS) = SMAX
         VV(IS) = VMAX
C  spline fit this segment (through IS); IOP(2) = 3 forces first
C     derivative to SCR(IS)=0 at right edge of fit.
         IOP(2) = 3
         SCR(IS) = 0.0D0
C  if less than 4 points in spline fit can't do 4-pt.  Lagrance
C     interpolation at left edge, therefore force derviative to zero
C     there.
         IF (IS .LT. 4) IOP(1) = 3
         CALL SPL1D1(IS, SV, VV, SCR, IOP, 1, AV, BV, CV)
         CALL SPL1B1(IS, SV, VV, SCR, 1, AV, BV, CV, DV)
         NS = IS
C        WRITE(6, 6600) NS, SMAX, (SV(I), I=1, NS)
C6600    FORMAT(' VSPLIN,NS,SMAX=', I5, F15.7/' SV='/(1X, 8F15.7))
C        WRITE(6, 6601) (VV(I), I=1,NS)
C6601    FORMAT(' VV='/(1X, 1P8E15.7))
C        WRITE(6, 6602) AV(1), BV(1), CV(1), DV(1)
C6602    FORMAT(' COEFF', 2X, 1P4E15.7)
C  set constants need for spline fit of next segment; derivative at
C     left egde of next segment will be force to zero.
         IORG = IS
         IOP(1) = 3
         SCR(1) = 0.0D0
         IS = ISS - 1
      END IF
C
      IF (NS .LT. NSLAST) THEN
C  more points left, spline fit next segment
         IF (LSYM .AND. NSLAST - IS .LT. 2) NSLAST = NSLAST + 2
C  load the rest of the SV and VV array; they are shifted from SS and
C     VAD by 1 if the maximum was added (IS0=1)
         I0 = IS + 1
         DO 30 I = I0,NSLAST
            IS = IS + 1
            ISP = IS + IS0
            SV(ISP) = SS(IS)
            VV(ISP) = VAD(IS)
   30    CONTINUE
         NN = ISP - IORG + 1
         IF (LSYM .AND. NSLAST.GT.NMID) THEN
            NS = NMID+IS0
         ELSE
            NS = ISP
         END IF
         IF (NN .LT. 4 .OR. (LSYM .AND. NSLAST.EQ.NMID)) THEN
            IOP(2) = 3
            SCR(NN) = 0.D0                                              GCL1092
         ELSE
            IOP(2) = 5
         END IF
         CALL SPL1D1(NN, SV(IORG), VV(IORG), SCR, IOP, 1, AV(IORG),
     *      BV(IORG), CV(IORG))
         CALL SPL1B1(NN, SV(IORG), VV(IORG), SCR, 1, AV(IORG), BV(IORG),
     *      CV(IORG), DV(IORG))
C        WRITE(6, 6603) NN, AV(1), BV(1), CV(1), DV(1)
C6603    FORMAT(' NN=', I5, /' COEFF', 2X, 1P4E15.7)
      END IF
      RETURN
      END
