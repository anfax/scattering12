!!!*************************************************************
! 文件/File: boltz.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: boltz.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE BOLTZ (BET, VM, PE, SUM, IFLG, ESV, EMAX, IEMAX, VMAX)
C
C     BOLTZ  - do numerical integration of P's to get kappas
C
C  Called by:
C     KAPVA  - compute kappas
C
C  Calls:
C     QUADFT - quadratic fit of three points
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION PE(2,1), SUM(4), ESV(1), EMAX(2), IEMAX(2)
      DIMENSION PEXP(2), PEXPMX(2), E(3), F(3), COEF(3)
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /QUADKA/ PT(NQGKDM), WT(NQGKDM,2), NQ12, NSEG
C
      DO 200 I = 1, 3                                                   GCL1096
             E(I) = 0.D0                                                GCL1096
200   CONTINUE                                                          GCL1096
C
      XN = 0.5D0*BET*VM
      XINSEG = 1.D0/NSEG
      IF (IFLG .NE. 0) XN = -XN
      DO 10 J = 1,4
         SUM(J) = 0.D0                                                  GCL1092
   10 CONTINUE
      PEXPMX(1) = 0.D0                                                  GCL1092
      PEXPMX(2) = 0.D0                                                  GCL1092
      NPE = 0
C  Loop over number of segments
      DO 60 L = 1,NSEG
         CL = 2.D0 - (2.D0*DBLE(L) - 1.D0)*XINSEG                       GCL1096
C  Loop over number of quadrature points in each segment
         DO 50 N = 1,NQ12
            NPE = NPE + 1
            T = EXP(XN*(CL-PT(N)*XINSEG))
C  Loop over the two Kronrod quadratures in S integration
            DO 20 K = 1,2
               P = PE(K,NPE)
               PEXP(K) = T*P
               IF (IFLG .NE. 0) P = 1.D0-P
               P = P * T
               IF (P .GE. PEXPMX(K)) THEN
                  PEXPMX(K) = P
                  IEMAX(K) = NPE
               END IF
   20       CONTINUE
C  Loop over the two Kronrod quadratures in the E integration
            DO 40 J = 1,2
               J0 = 2*(J-1)
               W = WT(N,J)
               DO 30 K = 1,2
                  II = K + J0
                  SUM(II) = SUM(II) + W*PEXP(K)
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      C0 = 0.D0                                                         GCL1092
      IF (IFLG .NE. 0) C0 = 1.D0                                        GCL1092
      DO 70 II = 1,4
         SUM(II) = C0+XN*SUM(II)*XINSEG
   70 CONTINUE
C  quadratic fit to maximum of P(E)*EXP(-BETA*E)
      DO 90 I = 1,2
         II = IEMAX(I)
         IF (II .EQ. 1) THEN
            IF (IFLG .EQ. 0) EMAX(I) = E(1)*CKCAL
            IF (IFLG .NE. 0) EMAX(I) = (2.0D0*VMAX-E(1))*CKCAL
         ELSE IF (II .EQ. NPE) THEN
            EMAX(I) = 0.0D0
         ELSE
            II = II - 1
            DO 80 K = 1,3
               T = PE(I,II+K-1)
               E(K) = ESV(II+K-1)/CKCAL
               IF (IFLG .NE. 0) THEN
                  T = 1.D0-T
                  E(K) = 2.D0*VMAX-E(K)
               END IF
               F(K) = T*EXP(-BET*E(K))
   80       CONTINUE
            CALL QUADFT(E, F, COEF)
            IF (COEF(3) .NE. 0.0D0) EMAX(I) = -0.5D0*CKCAL*COEF(2)/COEF(
     V3)
         END IF
   90 CONTINUE
      RETURN
      END
