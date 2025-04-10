!!!*************************************************************
! 文件/File: brnuli.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: brnuli.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE BRNULI(N,B)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION B(N)
C     B(I) CONTAINS THE 2I BERNOULI NUMBER,
C     I.E., B(1) IS B SUB 2 AS DEFINED IN ABRAM. AND STEGUN
      IF(N.GT.31) GO TO 1000
      TPI = 8.0D0*ATAN(1.0D0)
      ALTPI = LOG(TPI)
      CALL ALFCT(N,B)
      B(1) = 1.0D0/6.0D0
      B(2) = -1.0D0/30.0D0
      SGN = -1.0D0
      DO 100 I = 3,N
      SGN = -SGN
      J = 2*I
      FJ = DBLE(J)                                                      GCL1092
      T = 2.0D0*SGN*EXP(B(I) - FJ*ALTPI)
      SUM = 1.0D0
      DO 50 K = 2,400
      FK = DBLE(K)                                                      GCL1092
      TERM = FK**(-J)
      SUM = SUM + TERM
   50 IF(TERM/SUM .LT. 1.D-13) GO TO 60
      WRITE(6,600) K
   60 B(I) = T*SUM
  100 CONTINUE
      RETURN
 1000 WRITE(6,6000)
6000  FORMAT(/,1X,T5,'Error: 31 is the largest number of the Bernouli ',
     *               'numbers that can be ',
     *       /,1X,T12,'computed with this version of the code.')
      STOP 'BRNULI 1'                                                   GCL0992
600   FORMAT(/,1X,T5,'Warning: In BRNULI, the series did not converge ',
     *               'after ',I10,' terms')
      END
