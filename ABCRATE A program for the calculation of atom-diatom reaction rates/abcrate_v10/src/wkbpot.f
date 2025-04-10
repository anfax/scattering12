!!!*************************************************************
! 文件/File: wkbpot.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: wkbpot.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      FUNCTION WKBPOT (U, DVDU, IRANGE)
C
C     WKBPOT - potential for u motion for fixed s
C
C  Called by:
C     PHSINT - compute phase integrals needed for WKB quantization
C     VEXTR  - find extrema in potential in u coordinate (for WKB
C              quantization)
C     WKBTP  - find turning points in potential for u motion
C
C  Calls:
C     PEF    - evaluates potential
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /WKB1/   X0, Y0, UX, UY
C
      X = X0+UX*U
      Y = Y0+UY*U
      R2 = Y/CM2
      R1 = X - CM1*R2
      R3 = R1 + R2
c     write (99,*) u, r1, r2, r3
      IRANGE = 1
      VV = 1.D35
      DVDU = -1.D30
      IF (R1 .GE. 1.D-4 .AND. R2 .GE. 1.D-4 .AND. R3 .GE. 1.D-4) THEN
         CALL PEF (R1, R2, R3, V, D1, D2, D3, 1)                        GCL0893
         IF (V .LT. VV) THEN
            IRANGE = 0
            VV = V
            DX = D1 + D3
            DY = (D2+D3 - CM1*DX)/CM2
            DVDU = UX*DX + UY*DY
         END IF
      END IF
      WKBPOT = VV
      RETURN
      END
