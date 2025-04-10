!!!*************************************************************
! 文件/File: lcsa2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lcsa2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE LCSA2 (CSA, ALFMN, DS, DX, DY, UX, UY)
C
C     LCSA2  - compute cosine factors for LAG calculation
C
C  Called by:
C     LCSA   - compute cosine factors for LAG calculation
C
C  Calls:
C     LCSA3  - compute cosine factors for LAG calculation
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION CSA(2,4), ALFMN(2)
      DT2 = DX*DX + DY*DY
      TDV = UX*DY - UY*DX
      DO 10 I = 1,2
         CALL LCSA3 (CSA(I,1), ALFMN(I), DS, DT2, TDV)
   10 CONTINUE
      CALL LCSA3 (CSA(1,2), 0.0D0, DS, DT2, TDV)
      CALL LCSA3 (CSA(1,3), 0.5D0, DS, DT2, TDV)
      CALL LCSA3 (CSA(1,4), 1.0D0, DS, DT2, TDV)
      CSA(2,2) = CSA(1,2)
      CSA(2,3) = CSA(1,3)
      CSA(2,4) = CSA(1,4)
      RETURN
      END
