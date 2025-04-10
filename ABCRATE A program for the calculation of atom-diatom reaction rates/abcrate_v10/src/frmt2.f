!!!*************************************************************
! 文件/File: frmt2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: frmt2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE FRMT2 (FORM, IFORM, T)
C
C     FRMT2  - computes format for printing ratios of rates
C
C  Called by:
C     TRANS  - print out transmission coefficients (kappas)
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*200 FORM
      IF (T .LT. 95.D0) THEN
         IFF = IFORM + 9
         FORM(IFORM:IFF) = ',0PF5.0,2X'
      ELSE IF (T .LT. 99.5D0) THEN
         IFF = IFORM + 9
         FORM(IFORM:IFF) = ',0PF6.1,1X'
      ELSE
         IFF = IFORM + 6
         FORM(IFORM:IFF) = ',0PF7.2'
      END IF
      IFORM = IFF + 1
      RETURN
      END
