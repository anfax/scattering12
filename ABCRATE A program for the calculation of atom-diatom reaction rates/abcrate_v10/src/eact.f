!!!*************************************************************
! 文件/File: eact.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: eact.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE EACT (I0, N, NDIM, RATE, EAC, A, S, IDX)               TCA0197
C
C     EACT   - compute activation energy for temperature range T1 to T2
C
C  Called by:
C     SUMMRY - summarize rate constants
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NKATYP=21)                                             TCA0197
      PARAMETER (NOB=4)                                                 TCA0197
      PARAMETER (NKALAB=NKATYP-NOB)                                     TCA0197
      DIMENSION RATE(NDIM,N), EAC(N), A(N), S(2,N), IDX(NKALAB)         TCA0197
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /EACT1/   IACT, NT1(10), NT2(10)
      PARAMETER (NTEMDM=100)
      DOUBLE PRECISION KAPW
      LOGICAL PTEMP 
      COMMON /TEMPCM/ TEMP(NTEMDM), BETA(NTEMDM), CPHI(NTEMDM),
     *CNST(NTEMDM), CPHIC(NTEMDM), CNSTC(NTEMDM), EQUIL(NTEMDM,2),
     *KAPW(NTEMDM,2), RATIO(NTEMDM,2), NTMAX, PTEMP(NTEMDM)
C
      IF (IACT .LE. 0) RETURN
C  Loop over number of temperature pairs
      DO 20 IT = 1,IACT
         N1 = NT1(IT)
         N2 = NT2(IT)
         N2 = MIN(N2,NTMAX)
         N1 = MIN(N1,N2-1)
         N1 = MAX(1,N1)
         IF (N1 .EQ. N2) GO TO 20
         T1 = TEMP(N1)
         T2 = TEMP(N2)
         B1 = BETA(N1)
         B2 = BETA(N2)
         C = CKCAL/(B1-B2)
C  Loop over methods
         DO 10 I = I0,N
            EAC(I) = 0.0D0
            A(I) = 0.0D0
            S(1,I) = 0.0D0
            S(2,I) = 0.0D0
            IF (RATE(N2,IDX(I)) .LE. 0.D0 .OR. RATE(N1,IDX(I)) .LE.     TCA0197
     *          0.D0) GO TO 10                                          TCA0197
            EAC(I) = C*LOG(RATE(N2,IDX(I))/RATE(N1,IDX(I)))             TCA0197
            A(I) = EXP(B1*EAC(I)/CKCAL + LOG(RATE(N1,IDX(I))))          TCA0197
            T = TPI*A(I)*TAU
            S(1,I) = CKCAL*(LOG(T*B1) - 2.0D0)/CKAU
            S(2,I) = CKCAL*(LOG(T*B2) - 2.0D0)/CKAU
   10    CONTINUE
C  Print out activation energy, Arrhenius parameter and entropy of
C     activation.
         WRITE (6, 600) T1, (EAC(I),I=I0,N)
         WRITE (6, 601) T2, (A(I),I=I0,N)
         WRITE (6, 602) T1, (S(1,I),I=I0,N)
         WRITE (6, 602) T2, (S(2,I),I=I0,N)
   20 CONTINUE
      RETURN
  600 FORMAT (3X, F9.3, '  EA', 1P,9E12.4)                              GCL0992
  601 FORMAT (' TO', F9.3, '   A', 1P,9E12.4)                           GCL0992
  602 FORMAT (' S0(T=', F9.3, ')', 1P,9E12.4)                           GCL0992
      END
