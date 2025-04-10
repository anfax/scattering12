!!!*************************************************************
! 文件/File: elpart.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: elpart.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE ELPART (AX, NE, EL, NDGEN)
C
C     ELPART - compute electronic partition functions
C
C  Called by:
C     DATAIN - read in data and set up constants
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*2 AX                                                    TCA1296
      LOGICAL LSET
      DIMENSION AX(3), NE(5), EL(5,5), NDGEN(5,5), Q(5)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      PARAMETER (NTEMDM=100)
      DOUBLE PRECISION KAPW
      LOGICAL PTEMP
      COMMON /TEMPCM/ TEMP(NTEMDM), BETA(NTEMDM), CPHI(NTEMDM),
     *CNST(NTEMDM), CPHIC(NTEMDM), CNSTC(NTEMDM), EQUIL(NTEMDM,2),
     *KAPW(NTEMDM,2), RATIO(NTEMDM,2), NTMAX, PTEMP(NTEMDM)
      COMMON /ELFCT/  ELFACT(NTEMDM,2)
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      DIMENSION EFACTF(NTEMDM), EFACTR(NTEMDM)
      EQUIVALENCE (EFACTF, ELFACT), (EFACTR, ELFACT(1,2))
C
      DO 10 IT = 1,NTMAX
         BETA(IT) = CKAU/TEMP(IT)
         EFACTF(IT) = 0.D0                                              GCL1092
         EFACTR(IT) = 0.D0                                              GCL1092
   10 CONTINUE
      LSET = .FALSE.
      DO 30 I = 1,5
         EZ = EL(1,I)
         EL(1,I) = 0.D0                                                 GCL1092
         IF (NE(I) .LE. 0) GO TO 30
         LSET = .TRUE.
         N = NE(I)
         IF (N .LE. 1) GO TO 30
         DO 20 J = 2,N
            EL(J,I) = (EL(J,I)-EZ)/CKCAL
   20 CONTINUE
   30 CONTINUE
      IF (.NOT.LSET) RETURN
      IF (LGS(1)) WRITE (6, 600)
      IF (LGS(1)) THEN
         WRITE (6, 601) AX(1), (NDGEN(I,1), EL(I,1), EL(I,1)*EAU, EL(I,
     *      1)*CKCAL, EL(I,1)*CCM, I=1,NE(1))
         WRITE (6, 602) AX(2), AX(3), (NDGEN(I,2), EL(I,2), EL(I,2)*EAU,
     *      EL(I,2)*CKCAL, EL(I,2)*CCM, I=1,NE(2))
         WRITE (6, 603) AX(1), AX(2), AX(3), (NDGEN(I,3), EL(I,3), EL(I,
     *      3)*EAU, EL(I,3)*CKCAL, EL(I,3)*CCM, I=1,NE(3))
         WRITE (6, 602) AX(1), AX(2), (NDGEN(I,4), EL(I,4), EL(I,4)*EAU,
     *      EL(I,4)*CKCAL, EL(I,4)*CCM, I=1,NE(4))
         WRITE (6, 601) AX(3), (NDGEN(I,5), EL(I,5), EL(I,5)*EAU, EL(I,
     *      5)*CKCAL, EL(I,5)*CCM, I=1,NE(5))
         WRITE (6, 605) ((AX(I), I=1,3), J=1,3)
      END IF
      DO 70 IT = 1,NTMAX
         BET = BETA(IT)
         DO 60 I = 1,5
            T = NDGEN(1,I)
            N = NE(I)
            IF (N .GT. 1) THEN
               DO 50 J = 2,N
                  T = T + NDGEN(J,I)*EXP(-BET*EL(J,I))
   50          CONTINUE
            END IF
            Q(I) = LOG(T)
   60    CONTINUE
         EFACTF(IT) = Q(3) - Q(1) - Q(2)
         EFACTR(IT) = Q(3) - Q(4) - Q(5)
         IF (LGS(1)) WRITE(6, 606) TEMP(IT), Q, EFACTF(IT), EFACTR(IT)
   70 CONTINUE
      RETURN
  600 FORMAT (//,1X, 16('*'), 2X, 'Electronic energies and partition',
     *   ' functions', 1X, 16('*')// 44X, 'Electronic energy'/ 4X,
     *   'Species', 5X, 'Degen.', 8X, 'Hartree', 8X, 'eV', 10X,
     *   'kcal', 8x, 'cm**-1')
  601 FORMAT (/ 6X, A2, 9X, I2, 7X, 1P,4E13.5/ (17X, I2, 7X, 4E13.5))
  602 FORMAT (/ 6X, 2A2, 7X, I2, 7X, 1P,4E13.5/ (17X, I2, 7X, 4E13.5))
  603 FORMAT (/ 6X, 3A2, 5X, I2, 7X, 1P,4E13.5/ (17X, I2, 7X, 4E13.5)) 
  605 FORMAT (//, 35X, 'ln(partition functions)', 21X, 'ln(electronic',
     *   ' factors)'/ 4X, 'T(K)', 7X, A2, 10X, 2A2, 8X, 3A2, 9X, 2A2,
     *   9X, A2, 10X, 'Forward', 7X, 'Reverse')
  606 FORMAT(1X, F9.2, 1P,5E13.5, 2X, 2E13.5)     
      END
