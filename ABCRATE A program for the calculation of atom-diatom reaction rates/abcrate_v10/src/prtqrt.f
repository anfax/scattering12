!!!*************************************************************
! 文件/File: prtqrt.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: prtqrt.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PRTQRT (BETA, FB, AB, GB, QB, DELPHI, LP1)
C
C     PRTQRT - compute partition function of harm-quartic potential
C
C  Called by:
C     PFCNB  - compute bending part. fcn.
C
C  Calls:
C     EBEND  - compute bending energy levels
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LP1, LCONV
      LOGICAL LMAX
      COMMON /EBND1/   LMAX
      COMMON /STATE2/ DGBND, LSBEND, NBND1, NBND2
      SAVE EPS                                                          TCA1097
      DATA EPS /1.D-6/
C
      EZ = EBEND (0, FB, AB, GB)
      IF (LMAX) THEN
         QB = 80.D0                                                     GCL1092
         DELPHI = 999.0D0
      ELSE
         EB = 1.0D0/BETA
         E = MAX(EB, EZ)
         T = -0.5D0*FB
         AB24 = AB/24.0D0
         DELPHI = 999.0D0
         T2 = T*T + 4.0D0*E*AB24
         IF (T2 .GE. 0.0D0) THEN
            T = 0.5D0*(T+SQRT(T2))/AB24
            IF (T .GE. 0.D0) DELPHI = 2.0D0*SQRT(T)
         END IF
         IF (LSBEND .EQ. 0) THEN
            N = 0
            SUM = 1.D0                                                  GCL1092
   10       CONTINUE
               N = N + 1
               T = EXP(-BETA*(EBN(N, FB, AB, GB)-EZ))
               IF (LMAX) GO TO 20
               SUM = SUM + T
               LCONV = T/SUM .LT. EPS
            IF (N .LT. 100 .AND. .NOT.LCONV) GO TO 10
   20       CONTINUE
            IF (LP1 .AND. .NOT. LCONV) WRITE (6, 600) N, E, SUM, T
            QB = -BETA*EZ + LOG(SUM)
            QB = QB + QB
         ELSE
            EB = EBN(NBND1, FB, AB, GB)
            IF (LMAX) THEN
               QB = 80.D0                                               GCL1092
            ELSE
               IF (NBND1 .EQ. NBND2) EB = 2.D0*EB
               IF (NBND1 .NE. NBND2) EB = EB + EBN(NBND2, FB, AB, GB)
               IF (LMAX) THEN
                  QB = 80.D0                                            GCL1092
               ELSE
                  QB = -BETA*EB + LOG(DGBND)
               END IF
            END IF
         END IF
      END IF
      RETURN
  600 FORMAT (2X,T5,'Warning: In PRTQRT, N=', I5, ', E=', 1PE13.5, 
     *              ', SUM=', E13.5, ', LAST TERM=', E13.5)
      END
