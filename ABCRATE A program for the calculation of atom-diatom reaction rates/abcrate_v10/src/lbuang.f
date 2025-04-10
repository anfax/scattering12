!!!*************************************************************
! 文件/File: lbuang.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lbuang.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LBUANG(ISINT,S,XAB,YAB,X,Y,UX,UY,UALF,ADIF)
C    
C     LBUANG  - compute UALF and angle between normalized gradient
C               vector at s and vector from xmep at s to point
C               on the tunneling path
C
C  Called by:
C     LBETAS - solve for s of beta
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C   interpolate x,y,ux,uy at s
      CALL AITKN2(ISINT, S, X, Y, UX, UY)
C   compute UALF
      DELX = XAB - X
      DELY = YAB - Y
      UALF = SQRT (DELX*DELX + DELY*DELY)
C   compute angle
      ADIF = -UY*DELX + UX*DELY
      RETURN
      END
