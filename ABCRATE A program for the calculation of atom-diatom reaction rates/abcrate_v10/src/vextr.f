!!!*************************************************************
! 文件/File: vextr.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vextr.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

      SUBROUTINE VEXTR
C
C     VEXTR  - find extrema in potential in u coordinate (for WKB
C
C  Called by:
C     WKB    - compute WKB energy levels for stretch
C
C  Calls:
C     WKBPOT - potential for u motion for fixed s
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      COMMON /MASSES/ XMA, XMB, XMC, XM, XMP, XMU, XMUP, CM1, CM2, CM1P,
     *CM2P
      COMMON /MORBC/  DBC, XKBC, AMBC
      COMMON /WKB1/   X0, Y0, UX, UY
      COMMON /WKB2/   UMIN, VMIN, UMAX, VMAX
C
      UMIN = 0.D0                                                       GCL1092
      VMIN = WKBPOT (UMIN, DVDU, IRANGE)
      IF (IRANGE .NE. 0) THEN
         WRITE (6, 6001) X0, Y0, UX, UY, UMIN
         STOP 'VEXTR 1'
      END IF
      DU = 0.001D0
      U = UMIN - DU
      V1 = WKBPOT (U, DVDU, IRANGE)
      IF (IRANGE .NE. 0) THEN
         WRITE (6, 6001) X0, Y0, UX, UY, U
         STOP 'VEXTR 2'
      END IF
      IF (V1 .LE. VMIN) THEN
         WRITE(6, 6000) X0, Y0, UX, UY, U, VMIN, V1
         STOP 'VEXTR 3'
      END IF
      U = UMIN + DU
      V2 = WKBPOT (U, DVDU, IRANGE)
      IF (IRANGE .NE. 0) THEN
         WRITE (6, 6001) X0,Y0, UX, UY, U
         STOP 'VEXTR 4'
      END IF
      IF (V2 .LE. VMIN) THEN
         WRITE(6, 6000) X0, Y0, UX, UY, U, VMIN, V2
         STOP 'VEXTR 5' 
      END IF
   10 CONTINUE
         V1 = V2
         DU = 2.D0*DU
         U = UMIN + DU
         V2 = WKBPOT (U, DVDU, IRANGE)
      IF (V2 .GT. V1 .AND. V2 .LE. DBC .AND. IRANGE .EQ. 0 .AND.
     *   DU .LE. 1000.D0) GO TO 10
      U = UMIN
      V2 = VMIN
      DU = 0.1D0*DU
   20 CONTINUE
      IF (ABS(DU) .LE. 1.D-6) GO TO 40
   30    CONTINUE
            V1 = V2
            U = U + DU
            IF (U .LT. UMIN) THEN
               DU = 0.5D0*ABS(DU)
               U = UMIN + DU
            END IF
            V2 = WKBPOT (U, DVDU, IRANGE)
            IF (IRANGE .NE. 0) THEN
               U = U - DU
               V2 = V1
               DU = 0.5D0*DU
               GO TO 20
            END IF
            T = (V2-V1)/DU
            IF (ABS(T) .LT. 1.D-5) GO TO 40
         IF (V2 .GT. V1) GO TO 30
         DU = -0.5D0*DU
         GO TO 20
   40 CONTINUE
C found maximum
      VMAX = V2
      UMAX = U
      RETURN
6000  FORMAT(2X,T5,'Error: IN VEXTR for X,Y,UX,UY=', 1P,4E13.5, 
     *      /,2X,T12,'U,VMIN,V(U)=',1P,3E13.5) 
6001  FORMAT(2X,T5,'Error: In VEXTR for X,Y,UX,UY=', 1P,4E13.5, 
     *     /,2X,T12,'and U=', E13.5,', the bond coordinates ',
     *              'are out of range')
      END
