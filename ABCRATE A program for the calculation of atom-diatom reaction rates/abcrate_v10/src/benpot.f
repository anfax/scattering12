!!!*************************************************************
! 文件/File: benpot.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: benpot.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE BENPOT (PHI, R1S, R2S, TR12)
C
C     BENPOT - evaluates bending potential
C
C  Called by:
C     BEND - compute bending potential parameters
C
C  Calls:
C     PEF    - potential subroutine
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /BENDOP/ THETA(2), A11, A12, A21, A22, IBOPT
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
C
      IF (PHI .LT. 0.0D0 .OR. PHI .GT. PI) THEN                         GCL1092
         WRITE (6, 6000) R2, R2, IBOPT, PHI, PI
         STOP 'BENPOT 1'
      END IF
      R3 = SQRT(R1S+R2S-TR12*COS(PHI))
      CALL PEF (R1, R2, R3, V, D1, D2, D3, 1)                           GCL0893
      RETURN
6000  FORMAT(/,2X,T5,'Error: In BEND PHI is out of range: R1, R2 = ',
     *      2(1PE13.5,1X),/,2X,T12,'IBOPT, PHI, PI = ',I5,2(1PE13.5,1X))
      END
