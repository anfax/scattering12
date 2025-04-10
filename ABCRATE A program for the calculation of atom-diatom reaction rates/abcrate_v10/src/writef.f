!!!*************************************************************
! 文件/File: writef.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: writef.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE WRITEF (ITYP, IUNIT, FORM, ATYP)
C
C     WRITEF - write with FORM= construct
C
C  Called by:
C     TABL20 - print out table of ratios
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*10 ATYP
      CHARACTER*102 FORM
      IF (ITYP.EQ.1) WRITE (IUNIT, FMT=FORM)
      IF (ITYP.EQ.2) WRITE (IUNIT, FMT=FORM) ATYP
      RETURN
      END
