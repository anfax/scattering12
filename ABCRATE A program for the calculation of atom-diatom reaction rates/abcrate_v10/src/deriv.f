!!!*************************************************************
! 文件/File: deriv.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: deriv.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE DERIV (X, Y, DX, DY)
C
C     DERIV - compute first derivatives of potential w.r.t.  x,y
C        coordinates
C
C  Called by:
C     D2DX2  - compute matrix of second derivatives in Jacobi coor.
C     D2VDU2 - compute second derivative along u coordinate
C     GRAD   - follow gradient
C     RPATH  - follow reaction path and store MEP info on grid
C        (external to ROOT2D)
C     SADDLE - find saddle point(s) and do normal mode analysis
C        (external to ROOT2D)
C     STEP2  - follow gradient by choosing step to minimize potential
C     VMIN   - find minimum energy along u coordinate
C
C  Calls:
C     POT
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
C
      R2 = Y/CM2
      R1 = X - CM1*R2
      R3 = R1+R2
      CALL PEF (R1, R2, R3, V, D1, D2, D3, 1)                           GCL0893
      DX = D1 + D3
      DY = (D2+D3 - CM1*DX)/CM2
      RETURN
      END
