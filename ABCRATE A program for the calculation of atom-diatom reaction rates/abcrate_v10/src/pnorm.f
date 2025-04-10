!!!*************************************************************
! 文件/File: pnorm.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: pnorm.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PNORM (NEMAX,NQKPDM,ESV,EVAG,PE,PNRM,GAMSV,THET1)      TCA0997
C
C     PNORM  - normalizae LAG probabilities
C
C  Called by:
C     LAG    - compute LAG probabilities
C
C  Calls:
C     LPVAG  - extrapolate theta to E=VA
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LSAME1
      DIMENSION ESV(NQKPDM),PE(2,NQKPDM,8),PNRM(2,NQKPDM,8),
     *   GAMSV(2,NQKPDM),THET1(2,NQKPDM)
      DIMENSION THT(3)
C  normalize probabilities
C  RLCG
      DO 30 I = 1,2
         DO 10 IE = 1,3
            THT(IE) = -0.5D0 * LOG(PE(I,NEMAX-3+IE,1))
   10    CONTINUE
         CALL LPVAG (I, 1, EVAG, ESV(NEMAX-2), THT, TVAG)
         T1 = EXP(2.D0*TVAG)
         DO 20 IE = 1,NEMAX
            P1 = PE(I,IE,1)
            IF (P1 .GT. 0.D0) THEN
               PV = P1
               T = (1.D0+0.5D0*(T1-1.0D0)*P1*T1)*PV/(1.D0+PV)
               IF (T .LT. 1.0D0) THEN
                  PNRM(I,IE,1) = T
               ELSE
                  PNRM(I,IE,1) = P1
                  J = 1
                  WRITE (6,6000) I, J, IE, T
               END IF
            END IF
   20    CONTINUE
   30 CONTINUE
C  RLAG
      DO 80 I = 1,2
         DO 60 IE = 1,3
            THT(IE) = -0.5D0 * LOG(PE(I,NEMAX-3+IE,2))
   60    CONTINUE
         CALL LPVAG (I, 2, EVAG, ESV(NEMAX-2), THT, TVAG)
         T1 = EXP(2.D0*TVAG)
         DO 70 IE = 1,NEMAX
            P1 = PE(I,IE,2)
            IF (P1 .GT. 0.D0) THEN
               PV = P1
               T = (1.D0+0.5D0*(T1-1.0D0)*P1*T1)*PV/(1.D0+PV)
               IF (T .LT. 1.0D0) THEN
                  PNRM(I,IE,2) = T
               ELSE
                  PNRM(I,IE,2) = P1
                  J = 2
                  WRITE (6,6000) I, J, IE, T
               END IF
            END IF
   70    CONTINUE
   80 CONTINUE
C  LAG, LAG0.0, LAG0.5, LCG3
      DO 120 J = 1,4
         JJ = J + 2
         DO 110 I = 1,2
            DO 90 IE = 1,3
               THT(IE) = -0.5D0 * LOG(PE(I,NEMAX-3+IE,JJ))
   90       CONTINUE
            CALL LPVAG (I, JJ, EVAG, ESV(NEMAX-2), THT, TVAG)
            T1 = EXP(2.D0*TVAG)
            DO 100 IE = 1,NEMAX
               P1 = PE(I,IE,JJ)
               IF (P1 .GT. 0.D0) THEN
                  PV = P1
                  T = (1.D0+0.5D0*(T1-1.0D0)*P1*T1)*PV/(1.D0+PV)
                  IF (T .LT. 1.0D0) THEN
                     PNRM(I,IE,JJ) = T
                  ELSE
                     PNRM(I,IE,JJ) = P1
                     WRITE (6,6000) I, JJ, IE, T
                  END IF
               END IF
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
      RETURN
6000  FORMAT(2X,T5,'For I, J, IE=',3I4,', the normalized probability ',
     *             'was',1PE13.5,', but has been reset to 1.0')
      END
