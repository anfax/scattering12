!!!*************************************************************
! 文件/File: store.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: store.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE STORE (IS)
C
C     STORE  - store reaction path info
C
C  Called by:
C     ADIAB  - compute adiabatic potential curves
C     GTST   - compute free energies, CVT and ICVT rates
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      SS(IS) = S
      VS(IS) = VNOW
      DS(IS) = D
      XKS(IS) = XKM
      X2(IS)= X
      Y2(IS) = Y
      UXS(IS) = UX
      UYS(IS) = UY
      FBS(IS) = FB
      QFBS(IS) = QFB
      GBS(IS) = GB
      XMOMS(IS) = XMOM
      CAPS(IS) = CUR
      RETURN
      END
