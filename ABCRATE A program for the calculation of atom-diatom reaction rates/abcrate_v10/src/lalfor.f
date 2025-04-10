!!!*************************************************************
! 文件/File: lalfor.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lalfor.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LALFOR (ALF,THET,ALPH,THT,IALFMN,THETMN)
C
C     LALFOR - sort ALF, THET into ALPH, THT arrays.  ALPH is in 
C              increasing order.  Alpha at the minimum theta and the 
C              two closest alphas are kept along with their associated 
C              theta values.  
C
C  Called by:
C     LALFMN - find minimum theta(alf)
C
C     IALFMN = index of minimum theta in the THT array.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ALPH(3),THT(3),THET(2),THETMN(2)
C
      IF (THET(2).LE.THT(IALFMN)) THEN
         II = 2
         IALFMN = 2
         IF (ALF.LT.ALPH(2)) THEN
            ALPH(3) = ALPH(2)
            THT(3) = THT(2)
         ELSE
            ALPH(1) = ALPH(2)
            THT(1) = THT(2)
         END IF
         THETMN(1) = THET(1)
      ELSE IF (IALFMN.EQ.1) THEN
         IF (ALF.LT.ALPH(2)) THEN
            II = 2
            ALPH(3) = ALPH(2)
            THT(3) = THT(2)
         ELSE
            II = 3
         END IF
      ELSE IF (IALFMN.EQ.2) THEN
         IF (ALF.LT.ALPH(2)) THEN
            II = 1
         ELSE
            II = 3
         END IF
      ELSE
         IF (ALF.LT.ALPH(2)) THEN
            II = 1
         ELSE
            II = 2
            ALPH(1) = ALPH(2)
            THT(1) = THT(2)
         END IF
      END IF
      ALPH(II) = ALF
      THT(II) = THET(2)
      RETURN
      END
