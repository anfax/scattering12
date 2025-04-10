!!!*************************************************************
! 文件/File: reflec.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: reflec.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE REFLEC (IS, ISAD)
C
C     REFLEC - reflect MEP info across symmetry point
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Calls:
C     SHIFT  - shift r.p. info in grid
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      NSHLF = NSSP(ISAD)
      WRITE (6, 616)
      NSMAX = MIN(2*IS-1, NSDM)
      NMID = MIN(IS, NSDM1/2)
C  NMID = (NSMAX+1)/2
      IF (IS .NE. NSHLF) THEN
C  if CAPS(IS)=0 set CAPS(IS) = CAPS(IS-1)
         IF (ABS(CAPS(IS)) .LT. 1.D-6) CAPS(IS) = CAPS(IS-1)
C  curvature at saddle point
         CAPS(NSHLF) = CAPS(NSHLF-1) - (CAPS(NSHLF+1)-CAPS(NSHLF-1))*
     *      SS(NSHLF-1)/(SS(NSHLF+1)-SS(NSHLF-1))
      END IF
      ISHIFT = NMID - IS
      IF (ISHIFT .NE. 0) THEN
         ISGN = 1
         CALL SHIFT (IS, ISHIFT, ISGN, 0.D0, ISAD, NSSP)
      END IF
      ISMX = NMID - 1
      NSMAX = 2*NMID - 1
      S0 = 2.D0*SS(NMID)
      I2 = NSMAX + 1
      T1 = CM1/SQRT(CM1*CM1+CM2*CM2)
      UY = SQRT(0.5D0*(1.0D0-T1))
      UX = SQRT(0.5D0*(1.0D0+T1))
      UXS(NMID) = UX
      UYS(NMID) = UY
      T1 = UX*UX - UY*UY
      T2 = 2.D0*UX*UY
      DO 10 IS = 1,ISMX
         I2 = I2 - 1
C  I2 = NSMAX + 1 - IS
         SS(I2) = -SS(IS) + S0
         VS(I2) = VS(IS)
         DS(I2) = DS(IS)
         XKS(I2) = XKS(IS)
         R2 = Y2(IS)/CM2
         R1 = X2(IS)-CM1*R2
         X2(I2) = R2+CM1*R1
         Y2(I2) = CM2*R1
         UXS(I2) = T1*UXS(IS) + T2*UYS(IS)
         UYS(I2) = T2*UXS(IS) - T1*UYS(IS)
         FBS(I2) = FBS(IS)
         QFBS(I2) = QFBS(IS)
         GBS(I2) = GBS(IS)
         XMOMS(I2) = XMOMS(IS)
         CAPS(I2) = CAPS(IS)
   10 CONTINUE
      RETURN
  616 FORMAT (/,1X,T5,'Note: Symmetric reaction, the parameters on ',
     *                'the product side of the',
     *        /,1X,T11,'symmetric stretch line are found ',
     *                 'by reflection',/)
      END
