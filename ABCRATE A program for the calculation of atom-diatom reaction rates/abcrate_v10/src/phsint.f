!!!*************************************************************
! 文件/File: phsint.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: phsint.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE PHSINT (E, X1, X2, THETA, DNDE, IERR)
C
C     PHSINT - compute phase integrals needed for WKB quantization
C
C  Called by:
C     WKB    - compute WKB energy levels for stretch
C
C  Calls:
C     WKBTP  - find turning points in potential for u motion
C     WKBPOT - potential for u motion for fixed s
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      PARAMETER (NQWKDM=81)
      COMMON /QUADWK/ PTWK(NQWKDM), WTWK(NQWKDM), NQWK
C
      CALL WKBTP (E, X1, X2, IERR)
      IF (IERR .EQ. 1) THEN
         THETA = 0.0D0
         DNDE = 0.0D0
         IERR = 0
      ELSE
         XSUM = 0.5D0*(X1+X2)
         XDIF = 0.5D0*(X2-X1)
         THETA = 0.0D0
         DNDE = 0.0D0
         DO 10 IX = 1,NQWK
            X = XSUM + PTWK(IX)*XDIF
            T = E - WKBPOT (X, DVDU, IRANGE)
            IF (T .LE. 0.0D0) GO TO 10
            T = SQRT(T)
            THETA = THETA + WTWK(IX)*T
            DNDE = DNDE + WTWK(IX)/T
   10    CONTINUE
         T1 = SQRT(2.0D0*XMU)*XDIF/PI
         T2 = 0.5D0*T1
         THETA = T1*THETA
         DNDE = T2*DNDE
      END IF
      RETURN
      END
