!!!*************************************************************
! 文件/File: pfcnr.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: pfcnr.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PFCNR (BETA, MOMENT, ISIGMA, QCL, QQM)
C
C     PFCNR  - compute rotational part. fcn.
C
C  Called by:
C     DATAIN - read in data and set up constants
C     PFCN   - compute GTS partition functions and free energy
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MOMENT
      LOGICAL LCONV
      SAVE EPS                                                          TCA1097
      DATA EPS /1.D-8/
C
      BETAB = 0.5D0*BETA/MOMENT
      QCL = 1.D0/BETAB
      QQM = 1.D0                                                        GCL1092
      J = 0
C TJP1 = 2*J + 1
C FJJ = J*(J+1)
      TJP1 = 1.D0                                                       GCL1092
      FJJ = 0.D0                                                        GCL1092
   10 CONTINUE
         J = J + 1
         FJJ = FJJ + TJP1 + 1.D0                                        GCL1092
         TJP1 = TJP1 + 2.D0                                             GCL1092
         T = TJP1*EXP(-BETAB*FJJ)
         QQM = QQM + T
         LCONV = T/QQM .LT. EPS
      IF (J .LT. 400 .AND. .NOT.LCONV) GO TO 10
C  if not converged add on classical approx for J>400
      IF (.NOT.LCONV) THEN
         T = QCL*EXP(-BETAB*160000.D0)
         QQM = QQM + T
      END IF
      QQM = QQM / DBLE(ISIGMA)                                          TCA1097
      QQM = LOG(QQM)
      QCL = QCL / DBLE(ISIGMA)                                          TCA1097
      QCL = LOG(QCL)
      RETURN
      END
