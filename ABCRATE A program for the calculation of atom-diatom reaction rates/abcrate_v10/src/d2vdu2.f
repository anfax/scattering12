!!!*************************************************************
! 文件/File: d2vdu2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: d2vdu2.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE D2VDU2 (XS, YS, UX, UY, D2V)
C
C     D2VDU2 - compute second derivative along u coordinate.  Second
C        derivative is found by taking numerical first derivatives of
C        analytic first derivatives.
C
C  Called by:
C     DATAIN - read in data and set up constants
C     GEOM   - find equilibrium geometry for reactants and products
C     PARAM  - compute parameters for bound motion along reaction path
C
C  Calls:
C     DERIV  - compute derivatives of potential w.r.t. x,y coordinates
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EPS = 1.D-5
      ERR = 1.D35
      D2OLD = 0.0D0
      DEL = 0.1D0
      DO 10 I = 1,5
         DEL = DEL*0.1D0
         X = XS + UX*DEL
         Y = YS + UY*DEL
         CALL DERIV(X, Y, DX, DY)
         DV1 = DX*UX+DY*UY
         X = XS- UX*DEL
         Y=  YS - UY*DEL
         CALL DERIV(X, Y, DX, DY)
         DV2 = DX*UX+DY*UY
         D2 = 0.5D0*(DV1 - DV2)/DEL
         T1 = ABS(1.0D0-D2OLD/D2)
         IF (T1 .LE. ERR) THEN
            D2V = D2
            ERR = T1
         END IF
         D2OLD = D2
   10 CONTINUE
      IF (ERR .GT. EPS) WRITE(6, 600) ERR
      RETURN
600   FORMAT(/,2X,T5,'Warning, In D2VDU2 the second derivative ',
     *               'error is ',1PE15.5)
      END
