!!!*************************************************************
! 文件/File: cobend.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: cobend.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:41
!*************************************************************

      SUBROUTINE COBEND (N, NOM, FB, QFB, GB, E)
C
C     COBEND   - compute semiclassical eigenvalue of centrifugal
C        oscillator
C
C  Called by:
C     ADIAB  - compute adiabatic potential curves
C     GTST   - compute free energies, CVT and ICVT rates
C     PRTQRT - compute partition function of harm-quartic potential
C     RPHSUM - summarize reaction path info
C     SADDLE - find saddle point(s) and do normal mode analysis
C     TABL21 - print out table of GTS info
C     THRCOR - compute threshold corrections for ICVT
C     VBEND  - compute bending energy levels
C     VSPLN2 - spline fit of adiabatic potential in product channel
C
C  Calls:
C     COBINT - compute phase integrals needed for WKB quantization
C     COBTP  - find turning points in centrifugal oscillator potential
C     COBVEX   - find extrema in centrifugal oscillator potential
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LMAX
      COMMON /EBND1/   LMAX
      COMMON /CONST/  PI, TPI, EAU, CKAU, CCM, CKCAL, BOHR, TAU
      LOGICAL LGS(10)
      COMMON /LOGIC/  LGS
      SAVE EPS                                                          TCA1097
      DATA EPS/1.D-8/
C
      AB = 0.0D0                                                        GCL1092
      IF (LGS(4)) AB = QFB
c      write (6,*) ' cobend called with n,nom,fb,ab,gb='
c      write (6,*) n,nom,fb,ab,gb
C  get the extrema for the total bend potential (including centrifugal
C     term)
      CALL COBVEX (NOM,FB,AB,GB,V0,VMAX)
c      write (6,*) ' v0,vmax=', v0,vmax
      IF (LMAX) RETURN
      TTT = 0.001D0                                                     GCL1092
      DEL = MIN(TTT,VMAX-V0)
      XN = DBLE(N) + 0.5D0                                              GCL1092
      IF (AB.LT.0.0D0) THEN                                             GCL1092
C  for AB<0, potential has a maximum, energy level are quasibound, and
C     above the maximum no energy levels exist.  Check if the desired
C     energy level lies below the maximum energy.
         E = VMAX - DEL
C  get turning points
         CALL COBTP (E,Q1,Q2,NOM,FB,AB,GB)
c         write (6,*) ' ab<0, e,q1,q2=', e,q1,q2
C  compute phase integral
         CALL COBINT (E,Q1,Q2,THETA,DNDE,NOM,FB,AB,GB)
c         write (6,*) ' theta,dnde=', theta,dnde
         IF (THETA.LT.XN) THEN
            LMAX = .TRUE.
            RETURN
         END IF
      END IF
      E = V0 + DEL
      EMIN = V0
      EMAX = VMAX
C
C  loop over iterations in newton-raphson root search, 
C     initial guess is small energy above minimum of well
C     EMIN and EMAX bracket root and if newton step is
C     outside range a new energy of (EMIN+EMAX)/2 is used.
      IC = 0
c      write (6,*) ' ic,e,q1,q2,theta,dnde,emin,emax='
   10 CONTINUE
C     get turning points
         CALL COBTP (E,Q1,Q2,NOM,FB,AB,GB)
C     compute phase integral
         CALL COBINT (E,Q1,Q2,THETA,DNDE,NOM,FB,AB,GB)
C         write (6,6600) ic,e,q1,q2,theta,dnde,emin,emax
C6600     format(1x,i3,1p7e12.4)
         E2 = E                                                         GCL0992
         DELF = THETA - XN
         IF (DELF .GT. 0.D0) EMAX = MIN(E, EMAX)
         IF (DELF .LT. 0.D0) EMIN = MAX(E, EMIN)
         IC = IC + 1
         E = E - DELF/DNDE
         IF (E .LE. EMIN .OR. E .GE. EMAX)
     *      E = 0.5D0*(EMIN+EMAX)
      IF (ABS(DELF) .GT. EPS .AND. IC.LT.50) GO TO 10
C
      IF (IC .GE. 50) WRITE (6, 601) IC, E2, THETA                      GCL0992
      RETURN
601   FORMAT(/,2X,T5,'Warning: In COBEND, after ',I3,' iterations ',
     *               'E = ',1PE13.5,' and F(E) = ',1PE13.5)
      END
