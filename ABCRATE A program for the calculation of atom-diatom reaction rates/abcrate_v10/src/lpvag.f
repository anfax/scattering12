!!!*************************************************************
! 文件/File: lpvag.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lpvag.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LPVAG (I, J, EVAG, E, THT, TVAG)
C
C     LPVAG  - extrapolate theta to E=VA
C
C  Called by:
C     PNORM  - normalizae LAG probabilities
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION E(3), THT(3), B(3)
C
      T = (THT(2)*(EVAG-E(3)) - THT(3)*(EVAG-E(2)))/(E(2)-E(3))
      CALL QUADFT (E, THT, B)
      T1 = B(1) + EVAG*(B(2)+EVAG*B(3))
C     WRITE(6, 610) I, J, T, T1
      IF(ABS(T-T1) .GT. 1.D-2) WRITE(6, 610) I, J, T, T1
      TVAG = T1
      RETURN
  610 FORMAT (' For I=', I2, ', J=',I2, ' linear and quadratic',
     *   ' extrapolation of THETA to E=VAG are',1P,2E13.5)              GCL0992
      END
