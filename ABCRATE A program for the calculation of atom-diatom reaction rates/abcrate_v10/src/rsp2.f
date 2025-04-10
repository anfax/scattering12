!!!*************************************************************
! 文件/File: rsp2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: rsp2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE RSP2(F,RT,IND,VECT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(3),RT(2),VECT(2,2)
C        EIGENVALUES
      SUM = (F(1) + F(3))*.5D0
      DIF = (F(1) - F(3))*.5D0
      B = F(2)
      BSQ = B*B
      Q = SQRT(DIF*DIF+BSQ)
      XLAMP = SUM + Q
      XLAMM = SUM - Q
      IF(XLAMM.GT.XLAMP) GO TO 20
      IM = 1
      IP = 2
      RT(1) = XLAMM
      RT(2) = XLAMP
      GO TO 30
   20 IM = 2
      IP = 1
      RT(1) = XLAMP
      RT(2) = XLAMM
   30 IF(IND.EQ.0) RETURN
C        EIGENVECTORS
      T1 = -DIF/Q
      T2 = .5D0*(T1+1.D0)
      T1 = .5D0*(T1-1.D0)
      V2IM = SQRT(-T1)
      V2IP = SQRT(T2)
      VECT(1,IM) = (DIF-Q)*V2IM/B
      VECT(2,IM) = V2IM
      VECT(1,IP) = (DIF+Q)*V2IP/B
      VECT(2,IP) = V2IP
      RETURN
       END
