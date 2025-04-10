!!!*************************************************************
! 文件/File: step1.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: step1.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE STEP1 (X, Y, DX, DY, STEP, V0, DELS, ERR)
C
C     STEP1  - take step along a vector then find minimum normal to
C              vector
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Calls:
C     VMIN   - find minimum energy along u coordinate
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      X0 = X
      Y0 = Y
      X = X - STEP*DX
      Y = Y - STEP*DY
      XSV = X
      YSV = Y
      CALL VMIN (X, Y, DX, DY, STEP, V0, .FALSE., IERR)
      IF (IERR .NE. 0) THEN
         WRITE (6, 6000)
         STOP 'STEP1 1'
      END IF
      T1 = X-XSV
      T2 = Y-YSV
      ERR = SQRT(T1*T1 + T2*T2)
      T1 = X - X0
      T2 = Y - Y0
      DELS = SQRT(T1*T1 + T2*T2)
      RETURN
6000  FORMAT(/,2X,T5,'Error: In STEP1 there is a problem with VMIN')
      END
