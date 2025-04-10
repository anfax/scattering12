!!!*************************************************************
! 文件/File: vmin.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vmin.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE VMIN (XS, YS, DX, DY, DEL, V0, LROT, IERR)
C
C     VMIN   - find minimum energy along u coordinate
C
C  Called by:
C     CURVE  - compute curvature by finite difference
C     DATAIN - read in data and set up constants
C     GEOM   - find equilibrium geometry for reactants and products
C     GRAD   - follow gradient
C     INTERP - interpolate r.p. info from grid
C     RPATH  - follow reaction path and store MEP info on grid
C     STEP1  - take step along a vector then find minimum normal to
C              vector
C
C  Calls:
C     DERIV  - compute derivatives of potential w.r.t. x,y coordinates
C     PEF    - evaluate potential
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LROT
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
      SAVE EPS, EPSV, SMALL                                             TCA1097
      DATA EPS/1.0D-10/, EPSV/1.0D-12/, SMALL/1.0D-30/
C
C  Find min along potential cut defined by DX,DY
C  If LROT=.TRUE. then let DXS,DYS change with each new estimate
C  of the location of the minimum
C
      IC  = 0
      IERR = 0
      STEP = DEL
      X = XS
      Y = YS
      IF (LROT) THEN
         CALL DERIV (X, Y, DX1, DY1)
         DMAG = SQRT(DX1*DX1+DY1*DY1)
         DX1 = DX1/DMAG
         DY1 = DY1/DMAG
         AC = DX*DX1 + DY*DY1
         IF (AC .GE. 0.9D0) THEN
            DX = DX1
            DY = DY1
         ELSE
            DX = DX + DX1
            DY = DY + DY1
            DMAG = SQRT(DX*DX+DY*DY)
            DX = DX/DMAG
            DY = DY/DMAG
         END IF
      ELSE
         R2 = Y/CM2
         R1 = X - CM1*R2
         R3 = R1+R2
         CALL PEF (R1, R2, R3, V, D1, D2, D3, 1)                        GCL0893
      END IF
      VZ = V
      UZ = 0.0D0
      XZ = X
      YZ = Y
C
C     IDBG = 99
C     WRITE (IDBG, 9900) XS, YS, DX, DY, DEL, VZ, LROT
C9900 FORMAT (/ ' ENTER VMIN; XS, YS, DX, DY, DEL, VZ, LROT=',/
C    *   5X, 1P6E19.11, 5X, L1)
C     WRITE (IDBG, 9901)
C9901 FORMAT ( 2X, 'IC', 4X, 'STEP', 8X,
C    *   'VZ-V', 14X, 'X', 18X, 'Y', 14X, 'DX', 11X, 'DY')
C
C  Find minimum in potential normal to  gradient (DX,DY)
C  The gradient is rotated as the point moves for LROT = .TRUE.
      IFLAG = 0
      SCALE = 0.1D0
   10 CONTINUE
         U = UZ + STEP
         IC = IC + 1
         IF (IC .GT. 100) THEN
            IERR = 1
            GO TO 20
         END IF
         X = XZ + DY*U
         Y = YZ - DX*U
         CALL DERIV (X, Y, DX1, DY1)
C        WRITE (IDBG, 9902) IC, STEP, VZ-V, X, Y, DX, DY
C9902 FORMAT (1X, I3, 1P2E12.4, 2E19.11, 2E13.5)
         IF (VZ .GT. V) THEN
            IF (SCALE .LT. 0.5D0) THEN
               STEP = STEP/SCALE
               SCALE = 0.5D0
               STEP = STEP*SCALE
               GO TO 10
            END IF
            IF (LROT) THEN
               DMAG = ABS(VZ - V)
               IF (ABS(VZ) .GT. SMALL) DMAG = ABS(DMAG/VZ)
               IF (DMAG .LT. EPSV) GO TO 20
               DMAG = SQRT(DX1*DX1+DY1*DY1)
               DX1= DX1/DMAG
               DY1= DY1/DMAG
               AC = DX*DX1 + DY*DY1
               IF (AC .GE. 0.9D0) THEN
                  DX = DX1
                  DY = DY1
               ELSE
                  DX = DX + DX1
                  DY = DY + DY1
                  DMAG = SQRT(DX*DX+DY*DY)
                  DX = DX/DMAG
                  DY = DY/DMAG
               END IF
               XZ = X
               YZ = Y
               U = 0.0D0
            END IF
            UZ = U
            VZ = V
            IFLAG = 1
         ELSE
            STEP = -STEP
            IFLAG = IFLAG + 1
            IF (IFLAG .NE. 1) THEN
               STEP = ABS(STEP*SCALE)
               IFLAG = 0
               IF (STEP .LE. EPS) THEN
                  XS = XZ + DY*UZ
                  YS = YZ - DX*UZ
                  GO TO 20
               END IF
            END IF
         END IF
         GO TO 10
   20 CONTINUE
      V0 = VZ
      RETURN
      END
