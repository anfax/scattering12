!!!*************************************************************
! 文件/File: surfit.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: surfit.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE SURFIT (DX, DY, SGN)
C
C     SURFIT - calls PARAM
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Calls:
C     PARAM  - compute parameters for bound motion along reaction path
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
C  normal coordinates
      UX = -SGN*DY
      UY = SGN*DX
      CALL PARAM
      RETURN
      END
