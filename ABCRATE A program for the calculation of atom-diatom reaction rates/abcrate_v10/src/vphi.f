!!!*************************************************************
! 文件/File: vphi.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vphi.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE VPHI (RR1, RR2, PHI, VV)
C
C     VPHI   - compute bending potential as function of angle
C
C  Called by:
C     EXTRAS - print out extra info about GTS
C
C  Calls:
C     PEF    - evaluate potential
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
      R1 = RR1
      R2 = RR2
      R3 = R1+R2
      CALL PEF (R1, R2, R3, V, D1, D2, D3, 0)                           GCL0893
      VV = V
      R3 = SQRT(R1*R1 + R2*(R2-2.0D0*R1*COS(PHI)))
      CALL PEF (R1, R2, R3, V, D1, D2, D3, 0)                           GCL0893
      VV = V - VV
      RETURN
      END
