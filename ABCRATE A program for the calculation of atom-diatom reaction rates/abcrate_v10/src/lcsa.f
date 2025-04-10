!!!*************************************************************
! 文件/File: lcsa.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: lcsa.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************


      SUBROUTINE LCSA (CSAL, CSAR, ALFMN, SL, SR, LRFLC, IOP)
C
C     LCSA   - compute cosine factors for LAG calculation
C
C  Called by:
C     LAG    - compute LAG probabilities
C
C  Calls:
C     AITKN2 - do Aitken interpolation of X, Y, UX, UY
C     LCSA2  - compute cosine factors for LAG calculation
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LRFLC
      DIMENSION  CSAL(2,4), CSAR(2,4), ALFMN(2)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      LOGICAL LSYM
      COMMON /SYM/    LSYM, NMID
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      IS1 = 1
      IS2 = NSMAX
      CALL AITKN2 (IS1, SL, XL, YL, UXL, UYL)
      IF (LSYM .AND. LRFLC) THEN
         T = CM2*XL - (1.D0+CM1)*YL
         DX = (CM1-1.D0)*T/CM2
         DY = T
         DELS = 2.D0*(SS(NMID) - SL)
      ELSE
         CALL AITKN2 (IS2, SR, XR, YR, UXR, UYR)
         DX = XR - XL
         DY = YR - YL
         DELS = SR - SL
      END IF
      IF (IOP .NE. 2) CALL LCSA2 (CSAL, ALFMN, DELS, DX, DY, UXL, UYL)
      IF (IOP .GT.  1 .AND.  .NOT.LSYM) CALL LCSA2 (CSAR, ALFMN, DELS,
     *   DX, DY, UXR, UYR)
      RETURN
      END
