!!!*************************************************************
! 文件/File: interp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: interp.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE INTERP (ISI)
C
C     INTERP - interpolate r.p.  info from grid to location S using
C        quadratic fits
C
C  Called by:
C     ADIAB  - compute adiabatic potential curves
C     GTST   - compute free energies, CVT and ICVT rates
C
C  Calls:
C     PARAM  - compute parameters for bound motion along reaction path
C     QUADFT - quadratic fit of three points
C     VMIN   - find minimum energy along u coordinate
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LROT, LFIRST
      DIMENSION B(3)
      INCLUDE 'abc.inc'                                                 TCA0197
      PARAMETER (NSDM1=NSDM+1, NSDM4=NSDM+4, NSDM10=NSDM+10)
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /PARAM1/  S, VNOW, D, XKM, FB, QFB, GB, XMOM, X, Y, UX, UY,
     *CUR
      COMMON /PEFCM/  R1, R2, R3, V, D1, D2, D3                         GCL0992
      PARAMETER (NSADDM=4)
      LOGICAL LRSTRT
      COMMON /RPTHCM/ DEL, DELSV, R1ASY, R2ASY, ACASY, EPSASY, NSMAX,
     *NSSP(NSADDM), LRSTRT
      COMMON /BLKPCM/ SS(NSDM10), VS(NSDM10), DS(NSDM10), XKS(NSDM10),
     *FBS(NSDM10), QFBS(NSDM10), GBS(NSDM10), XMOMS(NSDM10),
     *X2(NSDM10), Y2(NSDM10), UXS(NSDM10), UYS(NSDM10), CAPS(NSDM10)
C
      SAVE LFIRST                                                       TCA1097
      DATA LFIRST /.TRUE./
      IF (LFIRST) THEN
         DXASY1 = -1.D0                                                 GCL1092
         DYASY1 = 0.D0                                                  GCL1092
         T = SQRT(CM1*CM1 + CM2*CM2)
         DXASY2 = -CM1/T
         DYASY2 = -CM2/T
         LFIRST = .FALSE.
      END IF
      SGN = SIGN(1.D0, S)
C  interpolate X,Y to values at S
      CALL QUADFT(SS(ISI), X2(ISI), B)
      X = B(1) + S*(B(2)+S*B(3))
      CALL QUADFT(SS(ISI), Y2(ISI), B)
      Y = B(1) + S*(B(2)+S*B(3))
C  interpolate UX,UY to values at S
      DO 10 I = 1,3
         DX = SGN*UYS(ISI+I-1)
         DY = -SGN*UXS(ISI+I-1)
         IF (SGN .LT. 0.D0) AC = DX*DXASY1 + DY*DYASY1
         IF (SGN .GT. 0.D0) AC = DX*DXASY2 + DY*DYASY2
         IF (AC .GT. ACASY) GO TO 20
   10 CONTINUE
   20 CONTINUE
      IF (AC .LE. ACASY) THEN
         CALL QUADFT(SS(ISI), UXS(ISI), B)
         UX = B(1) + S*(B(2)+S*B(3))
         CALL QUADFT(SS(ISI), UYS(ISI), B)
         UY = B(1) + S*(B(2)+S*B(3))
         CALL QUADFT(SS(ISI), CAPS(ISI), B)
         CUR = B(1) + S*(B(2)+S*B(3))
         LROT = .TRUE.
C        WRITE (6,6600) X,Y,UX,UY,CUR
C6600    FORMAT (' INTERPOLATED VALUES, X,Y,UX,UY,CUR=',1P5E13.5)
      ELSE
         UX = UXS(ISI+1)
         UY = UYS(ISI+1)
         CUR = CAPS(ISI+1)
         LROT = .FALSE.
C        WRITE (6, 6601) X,Y,UX,UY,CUR
C6601    FORMAT (' NO INTERPOLATION; X,Y,UX,UY,CUR=',1P5E13.5)
      END IF
      XSV = X
      YSV = Y
      DX = SGN * UY
      DY = -SGN * UX
      DXSV = DX
      DYSV = DY
C  find minimum in u potential
      CALL VMIN(X, Y, DX, DY, .001D0, VNOW, LROT, IERR)
      IF (IERR .NE. 0) THEN
         WRITE (6, 6001)
         STOP 'INTERP 1'
      END IF
C  check if the X,Y values changed too much
      T1 = XSV-X
      T2 = YSV - Y
      T = SQRT(T1*T1+T2*T2)
      IF (T .GE. 1.D-3) THEN
         X = XSV
         Y = YSV
         DX = DXSV
         DY = DYSV
         WRITE (6, 6000)
      END IF
C  get potential parameters
      UX = -SGN*DY
      UY = SGN*DX
C     WRITE (6, 6602) X,Y,UX,UY
C6602 FORMAT (' CALL PARAM, X,Y,UX,UY=', 1P4E13.5)
      CALL PARAM
      RETURN
 6000 FORMAT (/,2X,T5,'Warning: In interpolation of X, Y new value is',
     *    ' very different from old value')
 6001 FORMAT (/,2X,T5,'Error: In INTERP, problem with VMIN')
      END
