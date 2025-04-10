!!!*************************************************************
! 文件/File: qmpart.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: qmpart.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE QMPART (BETA, IS, S, D, XK, QVIB, LP, LWKB)
C
C     QMPART - compute part. fcn. for stretch
C
C  Called by:
C     PFCNST - compute stretching part. fcn.
C
C  Calls:
C     WKBINT - interpolate WKB energy levels from grid
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LP, LWKB, LCONV
C
C  vibrational part. fcn.
C  zero of energy is bottom of vibrational well
      NMAX = INT((XK-1.D0)/2.D0) + 1
      TIXK = 2.D0/XK
      TK = .5D0*TIXK
      E = D*TK*(2.D0-TK)
      IF (LWKB) CALL WKBINT (IS, S, 0, E, UL, UG, PER)
      EZ = E
      SUM = 1.0D0
      I = 1
      LCONV = .FALSE.
   10 CONTINUE
      IF (I .GE. NMAX .OR. LCONV) GO TO 20
         I = I + 1
         TK = TIXK + TK
         E = D*TK*(2.D0-TK)
         T = EXP(-BETA*(E-EZ))
         SUM = SUM + T
         LCONV = T/SUM  .LT.  1.D-8
         GO TO 10
   20 CONTINUE
      QVIB = -BETA*EZ + LOG(SUM)
      IF (LP .AND.  .NOT.LCONV) WRITE(6, 600) NMAX, QVIB, SUM, T
      RETURN
  600 FORMAT (2X,T5,'Warning: stretching partition function not ',
     *              'converged',/,2X,T14,'NMAX=',I5,5X,'QVIB=',1PE15.5,
     *        5X, 'SUM=', E15.5, 5X, 'LAST TERM=', E13.5)
      END
