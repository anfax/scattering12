!!!*************************************************************
! 文件/File: lalfst.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lalfst.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LALFST (X0, Y0)
C
C     LALFST - sets up for calculation of alf=0 path for LAG
C     LALF0  - computes x and y as function of beta
C
C  Called by:
C     LAG    - compute LAG probabilities
C     LAGTH  - compute barrier LAG penetration integral (theta) for
C              a given alf
C     LBETAS - solve for s of beta
C     LIADIB - evaluate adiabatic integrand
C
C  Calls:
C     AITKN2 - do Aitken interpolation of X, Y, UX, UY
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NQGKDM=81, NQKPDM=4*NQGKDM)
      LOGICAL LLAG, LLAGRS
      COMMON /LAGCOM/ PT3(NQGKDM), WT3(NQGKDM,2), NQ32, NSEG3, IOPTAU,  TCA0197
     *                LLAG, LLAGRS                                      TCA0197
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /THLAG1/ XM1, YM1, X1, Y1, UXM1, UYM1, UX1, UY1, SL, SR,
     *XDIF, XSUM, YDIF, YSUM
      SAVE NFINAL                                                       TCA1097
      DATA NFINAL /0/
      SAVE SN, XT, ISINT                                                TCA0197
C
      SN = 0.5D0 * (SL + SR)                                            GCL1092
      XT = 0.5D0 * (SR - SL)                                            GCL1092
      ISINT = 1
      CALL AITKN2(ISINT, SL, XM1, YM1, UXM1, UYM1)
      ISINT = NSMAX
      CALL AITKN2(ISINT, SR, X1, Y1, UX1, UY1)
      XSUM = 0.5D0 * (X1 + XM1)                                         GCL1092
      XDIF = 0.5D0 * (X1 - XM1)                                         GCL1092
      YSUM = 0.5D0 * (Y1 + YM1)                                         GCL1092
      YDIF = 0.5D0 * (Y1 - YM1)                                         GCL1092
      ISINT = 1
      RETURN
C
      ENTRY LALF0 (BET, X0, Y0)
      S = SN + BET * XT
      CALL AITKN2 (ISINT, S, X0, Y0, UX, UY)
      RETURN
      END
