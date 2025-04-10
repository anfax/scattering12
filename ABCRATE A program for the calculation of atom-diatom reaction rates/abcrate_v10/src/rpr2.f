!!!*************************************************************
! 文件/File: rpr2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rpr2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RPR2 (NS, SS, VS, DS, XKS, FBS, QFBS, GBS, XMOMS, X2,
     *   Y2, UXS, UYS, CAPS)
C
C     RPR2   - read in reaction path info
C
C  Called by:
C     RPREAD - read in reaction path info
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SS(NS), VS(NS), DS(NS), XKS(NS), FBS(NS), QFBS(NS),
     *GBS(NS), XMOMS(NS), X2(NS), Y2(NS), UXS(NS), UYS(NS), CAPS(NS)
C
      READ(1,100) SS
      READ(1,100) VS
      READ(1,100) DS
      READ(1,100) XKS
      READ(1,100) FBS
      READ(1,100) QFBS
      READ(1,100) GBS
      READ(1,100) XMOMS
      READ(1,100) X2
      READ(1,100) Y2
      READ(1,100) UXS
      READ(1,100) UYS
      READ(1,100) CAPS
      RETURN
  100 FORMAT (/, (1X, E19.10, 3E20.10))
      END
