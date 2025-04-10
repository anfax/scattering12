!!!*************************************************************
! 文件/File: init.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: init.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE INIT (ISAD)
C
C     INIT   - initialization of arrays for following reaction path
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /BENDTS/ FBTS(NSADDM), QFBTS(NSADDM), GBTS(NSADDM),
     *XMOMTS(NSADDM)
      COMMON /MORSTS/ DTS(NSADDM), XKTS(NSADDM), AMTS(NSADDM),
     *OMTS(NSADDM), OMIMG(NSADDM)
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
      COMMON /SADDL1/ VSP(NSADDM), R1SP(NSADDM), R2SP(NSADDM),
     *XSP(NSADDM), YSP(NSADDM), SVECT(2,NSADDM), UVECT(2,NSADDM),
     *NSAD, NSADMX(2)
C
      VNOW = VSP(ISAD)
      D = DTS(ISAD)
      XKM = XKTS(ISAD)
      FB = FBTS(ISAD)
      QFB = QFBTS(ISAD)
      GB = GBTS(ISAD)
      XMOM = XMOMTS(ISAD)
      X = XSP(ISAD)
      Y = YSP(ISAD)
      UX = UVECT(1,ISAD)
      UY = UVECT(2,ISAD)
      CUR = 0.D0                                                        GCL1092
      RETURN
      END
