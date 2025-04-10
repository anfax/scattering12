!!!*************************************************************
! 文件/File: titles.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: titles.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE TITLES (IRW, IO, ITLE)
C
C     TITLES - print out 3-line title
C
C read/write title from/to device IO
C IRW <= -1, formatted read, do not set date and time
C IRW =   0, formatted read, set date and time
C IRW =>  1, formatted write
C
C  Called by:
C     ADIAB  - compute adiabatic potential curves
C     DATAIN - read in data and set up constants
C     GTST   - compute free energies, CVT, ICVT AND CUS rates
C     KAPVA  - compute kappas
C     LAG    - compute LAG probabilities
C     PREAD  - read in info needed for LAG calculation
C     PWRITE - write out info needed to restart LAG calculation
C     RPATH  - follow reaction path and store MEP info on grid
C     RPREAD - read in reaction path info
C     RPWRIT - write out reaction path info
C     SUMMRY - summarize rate constants
C     TABL20 - print out table of ratios
C     TABL21 - print out table of GTS info
C     TABL22 - print out table of rates
C
C  Calls:
C     DATTIM - get time and date
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 CTITLE(80,3,3)
      IT = MAX (1, ITLE)
      IT = MIN (3, IT)
      IF (IRW .LT. 0) THEN
         READ (IO, 500) ((CTITLE(I,J,IT),I=1,80),J=1,3)
      ELSE IF (IRW .EQ. 0) THEN
         CALL DATTIM(CTITLE(1,1,IT))
         READ (IO, 500) ((CTITLE(I,J,IT),I=1,80),J=2,3)
      ELSE IF (IRW .GE. 1) THEN
         WRITE (IO, 600) ((CTITLE(I,J,IT),I=1,80),J=1,3)
      END IF
      RETURN
  500 FORMAT(80A1)
  600 FORMAT(1X,80A1)
      END
