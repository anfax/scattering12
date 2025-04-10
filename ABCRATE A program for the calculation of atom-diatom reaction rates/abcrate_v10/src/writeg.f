!!!*************************************************************
! 文件/File: writeg.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: writeg.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE WRITEG (FORM, LERR, OUTPUT, N)
C
C     WRITEG - write with FORM= construct
C
C  Called by:
C     TRANS  - print out transmission coefficients (kappas)
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LERR
      CHARACTER*200 FORM
      DIMENSION OUTPUT(N)
      WRITE (6, FMT=FORM, ERR=100) OUTPUT
      LERR = .FALSE.
      RETURN
  100 LERR = .TRUE.
      RETURN
      END
