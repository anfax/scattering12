!!!*************************************************************
! 文件/File: theta3.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: theta3.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE THETA3 (XMU, A, B, E, TH1, DTH1, TH2, DTH2)
C
C     THETA3 - phase integral for semiclassical eigenvalues of double
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
      SAVE PI, TPI                                                      TCA1097
      DATA PI, TPI /3.14159265359D0, 6.28318530718D0/
      RTEB = SQRT(E*B)
      RTEB2 = 2.0D0*RTEB
      B3PI = B*3.0D0*PI
      T1 = A + RTEB2
      IF (RTEB2 .LE. A) THEN
C  E below barrier
         QS = 2.0D0*RTEB2/T1
         CALL ELLIP(QS, EL1, EL2)
         TH1 = SQRT(XMU*T1)*(A*(EL2-EL1)+RTEB2*EL1)/B3PI
         DTH1 = SQRT(XMU/T1)*EL1/PI
         T2 = A - RTEB2
         RS = T2/T1
         CALL ELLIP(RS, EL1, EL2)
         TH2 = -2.0D0*SQRT(XMU*T1)*(A*EL2-RTEB2*EL1)/B3PI
         DTH2 = SQRT(XMU/T1)*2.0D0*EL1/PI
      ELSE
C  E above barrier
         T2 = RTEB2 - A
         SS = 0.5D0*T1/RTEB2
         CALL ELLIP(SS, EL1, EL2)
         TH1 = SQRT(XMU*RTEB)*(T2*EL1+2.0D0*A*EL2)/B3PI
         DTH1 = SQRT(XMU/RTEB)*EL1/TPI
         TS = 0.5D0*T2/RTEB2
         CALL ELLIP(TS, EL1, EL2)
         TH2 = 2.0D0*SQRT(XMU*RTEB)*(T1*EL1-2.0D0*A*EL2)/B3PI
         DTH2 = SQRT(XMU/RTEB)*EL1/PI
      END IF
      RETURN
      END
