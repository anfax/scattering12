!!!*************************************************************
! 文件/File: ders.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: ders.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE DERS (IS, S, DUM, DEDS)
C
C     DERS   - derivatives of morse turning point and zero pt. energy
C              w.r.t. S
C
C  Called by:
C     ADIAB  - compute adiabatic potential curves
C     LALFST - set up zeroth order path for LAG
C     PSAG   - compute SAG-type probabilites
C
C  Calls:
C     LOCS   - locate position of s in grid
C     QUADFT - quadratic fit of three points
C     UTP    - find outer turning point in u motion
C     WKBINT - interpolation WKB eigenvalue from grid
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION EZ(3), UM(3), B(3)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /STATE/  TNP1, LSTATE, NSTATE
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL1096
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
C  Locate grid index near S
      CALL LOCS (IS, S)
      IS = MAX(1, IS-1)
      IS = MIN(NSMAX - 2, IS)
      ISI = IS - 1
C  Loop over three grid points to obtain turning point and zero-point
C     energy.
      DO 10 J = 1,3
         ISI = ISI + 1
         D = DS(ISI)
         XK = XKS(ISI)
         EZ(J) = 2.0D0*D/XK
         IF (LGS(3)) EZ(J) = EZ(J)*(1.D0-0.5D0/XK)
         IF (LGS(7)) THEN
            CALL WKBINT(ISI, SS(ISI), 0, EZ(J), UL, UM(J), PER)
            IF (LSTATE .NE. 0 .AND. NSTATE .NE. 0) CALL WKBINT(ISI,
     *         SS(ISI), NSTATE, E, UL, UM(J), PER)
         ELSE
            UM(J) = UTP(D, XK, XMU)
         END IF
   10 CONTINUE
C  quadratic fit of turning points
      CALL QUADFT(SS(IS), UM, B)
      DUM = B(2) + 2.0D0*B(3)*S
C  quadratic fit of zero-point energy
      CALL QUADFT(SS(IS), EZ, B)
      DEDS = B(2) + 2.0D0*B(3)*S
      RETURN
      END
