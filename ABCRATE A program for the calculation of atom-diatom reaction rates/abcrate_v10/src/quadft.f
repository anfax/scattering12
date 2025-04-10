!!!*************************************************************
! 文件/File: quadft.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: quadft.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE QUADFT (X, Y, Z)
C
C     QUADFT - quadratic fit of three points
C
C  Called by:
C     ADIAB  - compute adiabatic potential curves
C     BEND   - compute bending potential parameters
C     BOLTZ  - do numerical integration of P's to get kappas
C     DERS   - derivatives of morse turning point and zero pt. energy
C              w.r.t. s
C     GTST   - compute free energies, CVT and ICVT rates
C     INTERP - interpolate r.p. info from grid
C     LALFMN - find minimum theta(alf)
C     LPVAG  - extrapolate theta to E=VA
C     STEP2  - follow gradient by choosing step to minimize potential
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(3), Y(3), Z(3)
      T = X(3) - X(1)
      T1 = (Y(3)-Y(1))/T
      CP = (T1 - (Y(2)-Y(1))/(X(2)-X(1)))/(X(3)-X(2))
      BP = T1 - CP*T
      Z(3) = CP
      Z(2) = BP-2.D0*CP*X(1)
      Z(1) = Y(1) - X(1)*(BP - CP*X(1))
      RETURN
      END
