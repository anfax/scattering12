!!!*************************************************************
! 文件/File: aitknf.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: aitknf.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:40
!*************************************************************

      DOUBLE PRECISION FUNCTION AITKNF(Y,FX,X,N)
C      This function subprogram is compatible with the UCC function.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FX(N+1),X(N+1)                                          TCA0197
      NP1 = N + 1
      DO 50 I = 1,N
      IP = I + 1
      DO 20 J = IP,NP1
   20 FX(J) = (FX(I)*(X(J)-Y) - FX(J)*(X(I)-Y))/(X(J)-X(I))
   50 CONTINUE
      AITKNF = FX(NP1)
      RETURN
      END
