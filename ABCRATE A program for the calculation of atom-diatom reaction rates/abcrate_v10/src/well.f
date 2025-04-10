!!!*************************************************************
! 文件/File: well.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: well.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE WELL (XWL, YWL, UXWL, UYWL, IS)
C
C     WELL   - print out info about local well in reaction path
C
C  Called by:
C     RPATH  - follow reaction path and store MEP info on grid
C
C  Calls:
C     D2DX2  - compute matrix of second derivatives in Jacobi coor.
C     RSP2   - diagonalize 2x2 matrix
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(3), ROOT(2), U(2,2)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),  GCL96
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      WRITE (6, 602) XWL, YWL
      CALL D2DX2 (XWL, YWL, F, ERR)
      CALL RSP2 (F, ROOT, 1, U)
      WRITE (6, 600) ROOT, ((U(I,J),J=1,2), I=1,2)
      DO 10 I = 1,2
         IF (ROOT(I) .LE. 0.0D0) GO TO 10
         OM = SQRT(ROOT(I)/XMU)
         WRITE (6, 601) I, OM, OM*EAU, OM*CKCAL, OM*CCM
   10 CONTINUE
      T1 = U(1,1)*UXS(IS) + U(2,1)*UYS(IS)
      T2 = U(1,2)*UXS(IS) + U(2,2)*UYS(IS)
      IF (ABS(T1) .LE. ABS(T2)) THEN
         SGN = SIGN(1.D0,T2)
         UXWL = SGN*U(1,2)
         UYWL = SGN*U(2,2)
      ELSE
         SGN = SIGN(1.D0,T1)
         UXWL = SGN*U(1,1)
         UYWL = SGN*U(2,1)
      END IF
      RETURN
  600 FORMAT (1X,T5,'Normal mode analysis of the well',
     * /,1X,T5,'Eigenvalues', T20, 1P,2E13.5,
     * /,1X,T5,'Eigenvectors',T20, 2E13.5,/,T20, 2E13.5)
  601 FORMAT (1X,T5,'Harmonic frequency for mode', I2,/,1X,T5,
     *   1PE13.5, ' Hartree', 1PE13.5, ' eV', E13.5, ' kcal',
     *   E13.5, ' cm**-1')
  602 FORMAT (1X,T5,'X=', 1PE13.5, ', Y=', E13.5)
      END
