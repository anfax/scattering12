!!!*************************************************************
! 文件/File: abc.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: abc.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:40
!*************************************************************

C
         SUBROUTINE prepef
C
         IMPLICIT DOUBLE PRECISION (A-H, O-Z)
         COMMON /PT1CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT2CM/ NDER, NDUM(21)
C
         CALL PREPOT
C 
         RETURN
         END
C
         SUBROUTINE pef (R1, R2, R3, PE, DEDR1, DEDR2, DEDR3, IDER)
C
         IMPLICIT DOUBLE PRECISION (A-H, O-Z)
         COMMON /PT1CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT2CM/ NDER, NDUM(21)
C
         R(1) = R1
         R(2) = R2
         R(3) = R3
         NDER = IDER
C
         CALL POT
C
         PE    = ENERGY 
         DEDR1 = DEDR(1)
         DEDR2 = DEDR(2)
         DEDR3 = DEDR(3)
C 
         RETURN 
         END
