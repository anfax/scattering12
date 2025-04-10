!!!*************************************************************
! 文件/File: frmt1.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: frmt1.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE FRMT1 (FORM, IFORM, T)
C
C     FRMT1  - computes format for printing kappas
C
C  Called by:
C     TRANS  - print out transmission coefficients (kappas)
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*200 FORM
      IF (T .EQ. 0.D0) THEN
         IFF = IFORM + 12
         FORM(IFORM:IFF) = ',5X,0PF3.0,3X'
      ELSE IF (T .GT. 0.1D0 .AND. T .LT. 10.0D0) THEN                   GCL1092
         IFF = IFORM + 9
         FORM(IFORM:IFF) = ',4X,0PF7.4'
      ELSE IF (T .GE. 10.0D0 .AND. T .LT. 100.0D0) THEN                 GCL1092
         IFF = IFORM + 12
         FORM(IFORM:IFF) = ',3X,0PF7.3,1X'
      ELSE
         IFF = IFORM + 7
         FORM(IFORM:IFF) = ',1PE11.4'
      END IF
      IFORM = IFF + 1
      RETURN
      END
