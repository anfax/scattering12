!!!*************************************************************
! 文件/File: dbwl.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: dbwl.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE DBWL (XMU, A, B, V0, NQ, EM, EP, ITYP)
C
C     DBWL   - compute semiclassical eigenvalue of harm.-quart.
C        potential.
C
C  Called by:
C     EBEND  - compute bending energy levels
C
C  Calls:
C     DBWLRT - performs root search for semiclassical eigenvalue of
C        harm.-quart. potential.
C
C  For SGN>0 (ITYP=2), V0=0.  For SGN<0 (ITYP=3), V0 = 0.25*A*A/B
C     A,B > 0
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DELSV = 0.1D0*V0
      IF (DELSV .LE. 0.D0) DELSV = 0.01D0
      E2 = EM
      XN = DBLE(NQ) - 0.5D0                                             GCL1092
      EMIN = 0.D0                                                       GCL1092
      EMAX = 1.D30
      DEL = DELSV
      SGN = -1.0D0
C  Solve for two roots for one NQ
      CALL DBWLRT (E2, XN, A, B, XMU, ITYP, EMIN, EMAX, DEL, SGN)
      EM = E2
      EMIN = E2
      EMAX = 1.D30
      DEL = DELSV
      SGN = 1.0D0
      CALL DBWLRT (E2, XN, A, B, XMU, ITYP, EMIN, EMAX, DEL, SGN)
      EP = E2
      RETURN
      END
