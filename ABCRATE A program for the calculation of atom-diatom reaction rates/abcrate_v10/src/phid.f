!!!*************************************************************
! 文件/File: phid.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: phid.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PHID (EPS, PH, DPH)
C
C     PHID   - compute phi function needed for uniform semiclassical
C              quantization of double well potential
C
C  Called by:
C     DBWL   - compute semiclassical eigenvalue of harm.-quart.
C        potential.
C
C  Calls:
C     BRNULI - computes Bernoulli numbers
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(10)
      LOGICAL LFIRST, LCONV
      SAVE LFIRST, NTERM, SWTCH, PSIHLF                                 TCA1097
      DATA LFIRST /.TRUE./, NTERM /10/, SWTCH /2.0D0/
      DATA PSIHLF /-1.963510026021423D0/
      AEPS = ABS(EPS)
      PH = 0.0D0
      DPH = 1.D35
      IF (EPS .EQ. 0.0D0) RETURN
      IF (AEPS .GE. SWTCH) THEN
         IF (LFIRST) THEN
            CALL BRNULI (NTERM, C)
            SGN = -1.0D0
            DO 10 I = 1,NTERM
               FJ = DBLE(2*I)                                           GCL1092
               SGN = -SGN
               C(I) = SGN*C(I)*(1.0D0-2.0D0**(1-2*I))/(FJ*(FJ-1.0D0))
   10       CONTINUE
            LFIRST = .FALSE.
         END IF
         T = (1.D0/AEPS)/AEPS
         SUM1 = C(1)
         SUM2 = -C(1)
         T2N = 1.D0                                                     GCL1092
         I = 1
   20    CONTINUE
            I = I + 1
            T2N = T2N*T
            TERM1 = C(I)*T2N
            SUM1 = SUM1 + TERM1
            TERM2 = -TERM1*(2.D0*I-1.D0)
            SUM2 = SUM2 + TERM2
            LCONV = ABS(TERM1/SUM1) .LT.  1.D-4 .AND.  ABS(TERM2/SUM2)
     *         .LT.  1.D-3
         IF (I .LT. NTERM .AND. .NOT.LCONV) GO TO 20
         IF (.NOT.LCONV) WRITE (6, 600) SWTCH, NTERM, SUM1, TERM1, SUM2,
     *      TERM2
         PH = SUM1/EPS
         DPH = SUM2*T
      ELSE
         SUM1 = 0.D0                                                    GCL1092
         SUM2 = 0.D0                                                    GCL1092
         TE = 2.0D0*AEPS
         I = 0
   30    CONTINUE
            I = I + 1
            TIM = 2.D0*DBLE(I) - 1.D0                                   GCL1092
            T = TE/TIM
            TERM1 = T - ATAN(T)
            SUM1 = SUM1 + TERM1
            T = T*T
            TERM2 = 2.D0*T/((1.D0+T)*TIM)
            SUM2 = SUM2 + TERM2
            LCONV = ABS(TERM1/SUM1) .LT.  1.D-6 .AND.  ABS(TERM2/SUM2)
     *         .LT.  1.D-6
         IF (I .LT. 800 .AND. .NOT.LCONV) GO TO 30
         IF (.NOT.LCONV) WRITE (6, 601) SWTCH, I, SUM1, TERM1, SUM2,
     *      TERM2
         SUM1 = SUM1*SIGN(1.D0, EPS)
         T = -LOG(AEPS) + PSIHLF
         PH = EPS*(1.D0+T) + SUM1
         DPH = T + SUM2
      END IF
      RETURN
  600 FORMAT(/,2X,T5,'Warning: In PHID PHI is not converged, EPS is ',
     *               'greater than ', F3.1,
     *       /,2X,T14,'NTERM,SUM1,TERM1,SUM2,TERM2=',I3,1P,4E13.4)
  601 FORMAT(/,2X,T5,'Warning: In PHID PHI is not converged, EPS is ',
     *               'less than ', F3.1, 
     *       /,2X,T14,'NTERM,SUM1,TERM1,SUM2,TERM2=',I3,1P,4E13.4) 
      END
