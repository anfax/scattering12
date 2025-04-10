!!!*************************************************************
! 文件/File: rpw2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rpw2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RPW2 (NS,SS,VS,DS,XKS,FBS,QFBS,GBS,XMOMS,X2,
     *Y2,UXS,UYS,CAPS)
C
C     RPW2   - write out reaction path info
C
C  Called by:
C     RPWRIT - write out reaction path info
C
C  Calls:
C     nothing
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SS(NS),VS(NS),DS(NS),XKS(NS),FBS(NS),QFBS(NS),
     *GBS(NS),XMOMS(NS),X2(NS),Y2(NS),UXS(NS),UYS(NS),CAPS(NS)
C
      WRITE(1,100) SS
      WRITE(1,101) VS
      WRITE(1,102) DS
      WRITE(1,103) XKS
      WRITE(1,104) FBS
      WRITE(1,105) QFBS
      WRITE(1,106) GBS
      WRITE(1,107) XMOMS
      WRITE(1,108) X2
      WRITE(1,109) Y2
      WRITE(1,110) UXS
      WRITE(1,111) UYS
      WRITE(1,112) CAPS
      RETURN
  100 FORMAT (' SS=',/,(1X,1PE19.10,3E20.10))
  101 FORMAT (' VS=',/,(1X,1PE19.10,3E20.10))
  102 FORMAT (' DS=',/,(1X,1PE19.10,3E20.10))
  103 FORMAT (' XKS=',/,(1X,1PE19.10,3E20.10))
  104 FORMAT (' FBS=',/,(1X,1PE19.10,3E20.10))
  105 FORMAT (' QFBS=',/,(1X,1PE19.10,3E20.10))
  106 FORMAT (' GBS=',/,(1X,1PE19.10,3E20.10))
  107 FORMAT (' XMOMS=',/,(1X,1PE19.10,3E20.10))
  108 FORMAT (' X2=',/,(1X,1PE19.10,3E20.10))
  109 FORMAT (' Y2=',/,(1X,1PE19.10,3E20.10))
  110 FORMAT (' UXS=',/,(1X,1PE19.10,3E20.10))
  111 FORMAT (' UYS=',/,(1X,1PE19.10,3E20.10))
  112 FORMAT (' CAPS=',/,(1X,1PE19.10,3E20.10))
      END
