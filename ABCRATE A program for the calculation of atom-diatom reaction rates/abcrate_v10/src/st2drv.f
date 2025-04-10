!!!*************************************************************
! 文件/File: st2drv.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: st2drv.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE ST2DRV (DX, DY, XZ, YZ, XD, YD, STEP, PHI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DX = COS(PHI)
      DY = SIN(PHI)
      XD = XZ - STEP*DX
      YD = YZ - STEP*DY
      CALL DERIV(XD, YD, DX, DY)
      RETURN
      END
