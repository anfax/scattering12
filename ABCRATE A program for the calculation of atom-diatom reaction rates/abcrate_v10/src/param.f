!!!*************************************************************
! 文件/File: param.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: param.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PARAM
C
C     PARAM  - compute parameters for bound motion along reaction path
C
C  Called by:
C     INTERP - interpolate r.p. info from grid
C     SURFIT - calls PARAM
C
C  Calls:
C     BEND   - compute bending potential parameters
C     D2VDU2 - compute second derivative along u coordinate
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /MORBC/  DBC, XKBC, AMBC
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
C  morse parameters
      CALL D2VDU2 (X, Y, UX, UY, D2)
      D = DBC - VNOW
      T = SQRT(XMU/D2)
      XKM = 4.D0*D*T
C  bending force constants and moment of inertia
      CALL BEND (X, Y, FB, QFB, GB, XMOM)
      RETURN
      END
