!!!*************************************************************
! 文件/File: curve.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: curve.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE CURVE (X0, Y0, DX0, DY0, DELC, SGN, CUR)
C
C     CURVE  - compute reaction path curvature by finite difference
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Calls:
C     VMIN   - find minimum energy along u coordinate
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CUR = 0.D0                                                        GCL1092
      IF (DELC .LE. 0.D0) RETURN
      STEP = 0.5D0*DELC
      DX1 = DX0
      DY1 = DY0
      X1 = X0 + STEP*DX1
      Y1 = Y0 + STEP*DY1
      CALL VMIN(X1, Y1, DX1, DY1, DELC, V0, .TRUE., IERR)
      IF (IERR .NE. 0) THEN
         WRITE (6, 6000)
         RETURN
      END IF
      DX2 = DX0
      DY2 = DY0
      X2 = X0 - STEP*DX2
      Y2 = Y0 - STEP*DY2
      CALL VMIN(X2, Y2, DX2, DY2, DELC, V0, .TRUE., IERR)
      IF (IERR .NE. 0) THEN
         WRITE (6, 6000)
         RETURN
      END IF
      T1 = X2 -X1
      T2 = Y2 -Y1
      DELS = SGN*SQRT(T1*T1 + T2*T2)
      T1 = (DX2-DX1)/DELS
      T2 = (DY2-DY1)/DELS
      CUR = T1*DY0 - T2*DX0
      RETURN
6000  FORMAT(/,2X,T5,'Warning: In CURVE, there is an error with VMIN so'
     *       /,2X,T14,'the curvature has been set equal to zero')
      END
