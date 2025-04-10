!!!*************************************************************
! 文件/File: theta2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: theta2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE THETA2 (XMU, A, B, E, TH1, DTH1, TH2, DTH2)
C
C     THETA2 - phase integral for semiclassical eigenvalues of double
C              well
C
C  Called by:
C     DBWL   - compute semiclassical eigenvalue of harm.-quart.
C        potential.
C
C  Calls:
C     ELLIP  - elliptic integrals of first and second kind
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE PI                                                           TCA1097
      DATA PI /3.14159265359D0/
      AH = 0.5D0*A
      TS = AH*AH + E*B
      T = SQRT(TS)
      SS = 0.5D0*(1.0D0-AH/T)
      B3PI = B*3.D0*PI
      T1 = 2.D0*SQRT(XMU*T)/B3PI
      CALL ELLIP(SS, EL1, EL2)
      TS = 0.5D0*(1.D0+AH/T)
      CALL ELLIP(TS, EL3, EL4)
      TH1 = T1*((AH+T)*EL1 - A*EL2)
      TH2 = 2.D0*T1*((T-AH)*EL3 + A*EL4)
      T1 = SQRT(XMU/T)/PI
      DTH1 = 0.5D0*T1*EL1
      DTH2 = T1*EL3
      RETURN
      END
