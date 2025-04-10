!!!*************************************************************
! 文件/File: cobprt.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: cobprt.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE COBPRT (BETA, FB, AB, GB, QB, DELPHI)
C
C     COBPRT - compute partition function of centrifugal osc.
C
C  Called by:
C     PFCNB  - compute bending part. fcn.
C
C  Calls:
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LMAX
      COMMON /COBND/  NLEVEL
      COMMON /EBND1/  LMAX
      LOGICAL LGS2
      COMMON /LOGIC2/ LGS2(10)
      COMMON /STATE2/ DGBND, LSBEND, NBND1, NBND2
C
      IF (LSBEND .EQ. 0) THEN
C  not state selected, sum over Boltzmann factors
C     get ground state energy
         CALL COBEND (0,0,FB,AB,GB,EZ)
         IF (LMAX) THEN
            QB = 80.D0                                                  GCL1092
            DELPHI = 999.0D0
         ELSE
            SUM = 1.0D0                                                 GCL1092
            EBMX = EZ
C     Sum over n and K where n = 2v + K
            DO 10 N = 1,NLEVEL
               K = N+2
C     Loop over K
5              CONTINUE
                  K = K - 2
                  NV = (N-K)/2
C     compute centrifugal oscillator energy level for bend
                  CALL COBEND(NV, K, FB, AB, GB, EB)
                  IF (.NOT.LMAX) THEN
                     DEGEN = 2.0D0                                      GCL1092
                     IF (K.EQ.0) DEGEN = 1.0D0                          GCL1092
                     T = EXP(-BETA*(EB-EZ))
                     SUM = SUM + DEGEN*T
                     EBMX = MAX(EB,EBMX)
                  END IF
               IF (K.GT.1) GO TO 5
10          CONTINUE
            TEST = EXP(-BETA*(EBMX-EZ))
            IF (TEST/SUM .GT. 1.D-6) WRITE (6, 6000) TEST,SUM           GCL1092
            QB = -BETA*EZ + LOG(SUM)
         END IF
      ELSE
C  state selected
C     compute centrifugal oscillator energy level for bend
         CALL COBEND(NBND1, NBND2,FB, AB, GB,EB)
         IF (LMAX) THEN
            QB = 80.D0                                                  GCL1092
         ELSE
            QB = -BETA*EB + LOG(DGBND)
         END IF
      END IF
      RETURN
6000  FORMAT(/,2X,T5,'Warning: In COBPRT, the calculation of the ',
     *               'partition fcn did not converge',
     *       /,2X,T14,'The last term and the sum = ',2(1PE15.7,1X))
      END
