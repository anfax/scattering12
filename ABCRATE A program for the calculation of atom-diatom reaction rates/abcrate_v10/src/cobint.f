!!!*************************************************************
! 文件/File: cobint.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: cobint.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE COBINT (E,Q1,Q2,THETA,DNDE,NOM,FB,AB,GB)
C
C     COBINT - compute phase integrals needed for semiclassical
C              quantization of the centrifugal oscillator
C
C  Called by:
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      PARAMETER (NQCODM=81)
      COMMON /QUADCO/ PTCO(NQCODM), WTCO(NQCODM), NQCO
C
      QSUM = 0.5D0*(Q1+Q2)
      QDIF = 0.5D0*(Q2-Q1)
      THETA = 0.0D0
      DNDE = 0.0D0
      DO 10 IQ = 1,NQCO
         Q = QSUM + PTCO(IQ)*QDIF
         QSQ = Q*Q
         T = 2.0D0*(E - 0.5D0*QSQ*(FB + AB*QSQ/12.0D0))/GB              GCL1092
         IF (NOM.NE.0) T = T - DBLE(NOM*NOM)/QSQ                        TCA0497
         IF (T .LE. 0.0D0) GO TO 10
         T = SQRT(T)
         THETA = THETA + WTCO(IQ)*T
         DNDE = DNDE + WTCO(IQ)/T
10    CONTINUE
      T1 = QDIF/PI
      T2 = T1/GB
      THETA = T1*THETA
      DNDE = T2*DNDE
      RETURN
      END
