!!!*************************************************************
! 文件/File: alfct.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: alfct.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:40
!*************************************************************

      SUBROUTINE ALFCT(N,A)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION A(N)
      A(1) = LOG(2.0D0)
      IF(N.LT.2) RETURN
      DO 100 I = 2,N
      J = 2*I
      FJ = DBLE(J)                                                      GCL1092
  100 A(I) = LOG(FJ) + LOG(FJ-1.0D0) + A(I-1)
      RETURN
      END
