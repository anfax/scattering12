!!!*************************************************************
! 文件/File: lin2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lin2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LIN2 (A11, A12, A21, A22, B1, B2, X1, X2)
C
C     LIN2   - solve two simultaneous linear equations
C
C  Called by:
C     BEND   - compute bending potential parameters
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      BOT=A22-(A21*A12/A11)
      TOP=B2-(A21*B1/A11)
      X2=TOP/BOT
      TOP=B1-A12*X2
      X1=TOP/A11
      RETURN
      END
