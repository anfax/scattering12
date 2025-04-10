!!!*************************************************************
! 文件/File: adda.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: adda.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

subroutine coulomb(atom,l,E,R,ac,f,h,fp,hp)
!program coulomb

! ADQAFGH.  FGH, A CODE FOR THE CALCULATION OF COULOMB RADIAL WAVE FUNCTIONS
! 1   FROM SERIES EXPANSIONS.  M.J. SEATON.                                 
! REF. IN COMP. PHYS. COMMUN. 146 (2002) 250                               
!***** code fgh.f
!  CALCULATION OF COULOMB FUNCTIONS AND
!  DERIVATIVES FROM SERIES SOLUTIONS.
!
!  FOR ATOM.EQ..TRUE.
!    READ L,E,R,AC
!    GETS FUNCTIONS f, g  AND h, AND THEIR DERIVATIVES
!  ELSE
!    READ L,ETA,RHO,AC
!    GETS FUNCTIONS F_l AND G_l, AND THEIR DERIVATIVES
!
!  FOR EACH FUNCTION y, THE QUANTITY [y] GIVES AN INDICATION OF 
!  THE NUMBER OF DECIMAL PLACES LOST DUE TO CANCELLATION EFFECTS
!  
use nrtype, only : dbl,i4b
implicit real(kind=dbl) (a-h,o-z)
      real(kind=dbl),PARAMETER :: PI=3.141592653589793238462643d0
      INTEGER(kind=i4b),PARAMETER :: nx=100
      LOGICAL ATOM 
      real(kind=dbl) E,R,ac,f,fp,h,hp
      integer(kind=i4b) :: l
!  
!write(6,*)'dentro',atom,l,E,R,ac
!------------------------
! Subroutine for evaluation of coulomb functions for matching at boundary with
! coulomb potential (ionic core)
!------------------------

!*******************************
!    1   READ(nx,500,END=2)ATOM
!READ(nx,500,END=2)ATOM
1 continue
!
IF(ATOM)THEN
!   READ(nx,*)L,E,R,AC
ELSE
!   READ(nx,*)L,ETA,RHO,AC
eta=E
   E=1./ETA**2
r=-eta*r
!           R=-ETA*RHO
ENDIF
!*******************************
!
CALL FGH(E,L,R,AC,F,FP,G,GP,H,HP,FBAR,GBAR,HBAR)
!
!write(6,*)'dentro2'
!-------------------------
! Output on screen
!WRITE(6,665)L,E,R,AC
!WRITE(6,666)F,FBAR,FP,G,GBAR,GP,H,HBAR,HP
!-------------------------
!
IF(.NOT.ATOM)THEN
  ETA=1./SQRT(E)
  D=1.
  DO J=1,L
    D=D*(1.+E*J**2)
  ENDDO
  D=D/(EXP(2.*PI*ETA)-1.)
  C=PI/(2.*ETA)
  CF=SQRT(C*D)*(-1)**(L+1)
  CG=SQRT(C/D)*(-1)**L
  BF=CF*F
  BG=CG*H
  BFP=-ETA*CF*FP
  BGP=-ETA*CG*HP
!---------------------------
! Output on screen
!  WRITE(6,667)L,ETA,RHO,AC
!  WRITE(6,668)BF,FBAR,BFP,BG,HBAR,BGP  
!---------------------------
ENDIF
!
  500   FORMAT(L2)
  665   FORMAT(' L=',I2,', E=',1P,E15.8,', R=',E15.8,&
     &  ', AC=',E9.2)
  666   FORMAT(&
     &  ' f=',1P,E15.8,' [',0P,F4.1,']',', fp=',1P,E15.8&
     &  ' g=',1P,E15.8,' [',0P,F4.1,']',', gp=',1P,E15.8&
     &  ' h=',1P,E15.8,' [',0P,F4.1,']',', hp=',1P,E15.8)
  667   FORMAT(/' L=',I2,', ETA=',1P,E15.8,', RHO=',E15.8,&
     &  ', AC=',E9.2)
  668   FORMAT(' F =',1P,E15.8,' [',0P,F4.1,']',', FP=',1P,E15.8&
     &         ' G =',E15.8,' [',0P,F4.1,']',', GP=',1P,E15.8)
!
!GOTO 1
!    2   STOP
!write(6,*)'end subroutine adda'
END subroutine coulomb
!end program coulomb
!******************************
!
      SUBROUTINE FGH(E,L,R,AC,F,FP,G,GP,H,HP,FBAR,GBAR,HBAR)
use nrtype
implicit real(kind=dbl)(a-h,o-z)
!
!  COULOMB FUNCTIONS f, g and h FROM POWER SERIES EXPANSIONS.
!
!      PARAMETER(PI=3.141592653589793238462643)
!
!  GET f AND g AND THEIR DERIVATIVES.
!
      CALL COULFG(L,E,R,AC,F,FP,G,GP,K,IERR,ACT,ABSF,ABSG)
!
      IF(IERR.EQ.2)THEN
!-----------------------
! Output on screen
!        WRITE(6,600)L,E
!-----------------------
        STOP
      ENDIF
!-----------------------
! Output on screen
!      IF(IERR.EQ.1)WRITE(6,610)L,E,K,ACT
!-----------------------
!
!  GET THE COULOMB FUNCTION h
!
!    CASE OF E.GE.0
      IF(E.GE.0)THEN
!       CALCULATE CAP. A
        A=1.
        IF(L.GT.0)THEN
          A1=1.
          A2=-E
          A3=E+E
          DO 10 I=1,L
            A2=A2+A3
            A1=A1+A2
   10     A=A*A1
        ENDIF
!       CALCULATE SCRIPT G 
        SG=A*FKHI(E,L,AC)/PI
        H=-G-SG*F
HP=-GP-SG*FP
!
!   CASE OF E.LT.0
      ELSE
        CALL ABG(E,L,AC,A,BG)
        H=-(G+BG*F)
        HP=-(GP+BG*FP)
      ENDIF
!
!  OBTAIN [y] QUANTITIES
absh=absg+ABS(sg)*absf
hbar=absh/ABS(h)-1
fbar=absf/ABS(f)-1
gbar=absg/ABS(g)-1
IF(hbar.GT.1)THEN
   hbar=LOG10(hbar)
ELSE
   hbar=0.
ENDIF
IF(fbar.GT.1)THEN
  fbar=LOG10(fbar)
ELSE
  fbar=0.
ENDIF
IF(gbar.GT.1)THEN
  gbar=LOG10(gbar)
ELSE
  gbar=0.
ENDIF
!
  600 FORMAT(///10X,'SERIES IN COULFG NOT CONVERGED'/&
     & 10X,'SMALL L = ',I5,', EPS = ',1PE12.4//)
  610 FORMAT(//5X,'*** FUNCTIONS FROM COULFG INACCURATE ***'/&
     & 5X,'SMALL L=',I2,', EPS = ',1P,E11.4,', K = ',I4,&
     & ',ACTACC = ',E9.2//)
      RETURN
      END 
!****************************************************8
!
      SUBROUTINE COULFG(LL,EPS,RHO,ACC,F,FP,G,GP,K,IERR,ACTACC,&
     & absf,absg)
use nrtype,only : dbl,i4b,Pi
implicit real(kind=dbl)(a-h,o-z)
!
!  CALCULATES COULOMB FUNCTIONS F AND G AND THEIR DERIVATIVES
!
!  ORIGINAL VERSION PUBLISHED IN COMP. PHYS. COMM.25, 87, 1982.
!  PRESENT VERSION MODIFIED TO AVOID UNDERFLOW AND OVERFLOW
!  CONDITIONS IN THE SUMMATIONS OVER N OF
!         U(N)=A(N)*RHO**(N+L+1)
!     AND V(N)=D(N)*RHO**(N+L+1)
!  U(N) AND V(N) ARE CALCULATED RECURSIVELY.
!
!  INPUT -
!        LL=ANGULAR MOMENTUM QUANTUM NUMBER
!        EPS=Z-SCALED ENERGY IN RYDBERGS
!        RHO=Z-SCALED RADIAL VARIABLE IN ATOMIC UNITS
!        ACC=ACCURACY REQUIRED
!
!  OUTPUT -
!        F=REGULAR FUNCTION
!        FP=DERIVATIVE OF F
!        G=IRREGULAR FUNCTION
!        GP=DERIVATIVE OF G
!        K=NUMBER OF TERMS NEEDED IN EXPANSION
!        IERR=ERROR CODE
!        ACTACC=ACCURACY ACTUALLY ACHIEVED
!
!  CONVERGENCE CRITERION -
!        VALUE OF WRONSKIAN CONVERGED TO ACCURACY OF 0.5*ACC
!
!  ERROR CODES -
!        IERR=0, CONVERGED WITH ACTACC.LT.ACC
!        IERR=1, CONVERGED WITH ACTACC.GT.ACC
!        IERR=2, NOT CONVERGED WITH 101 TERMS IN MAIN SUMMATION
!
!  INITIALIZATION
!
     real(kind=dbl),PARAMETER ::GAMMA=-0.577215664901532860606512d0,ONE=1.0d0,TWO=2.0d0,R2PI=ONE/(2.d0*PI),PS0=ONE+TWO*GAMMA
     real(kind=dbl) ::sp
!
      IERR=0
      LP1=LL+1
      L2=2*LL
      L2P1=L2+1
      FL=LL
      FLP1=LP1
      FL2P1=L2P1
      E2=.5*EPS
      R2=2.*RHO
      ACC2=2.*ACC
!
!     INITIALIZE FA=FACTORIAL(2*LL+1)
!     AND PS=PSI(2*LL+2)+PSI(1)
!
      FA=1.
      PS=PS0
!
!
!  CALCULATE ALPHA(N) AND BETA(N) AND INITIALIZE S AND SP
!  CONTINUE CALCULATION OF FA AND PS
!
!     S AND SP FOR N=0
      X3=-L2
      X2=L2P1
      X1=-2.*R2**(-LP1)
      SP=X3*X1
      X1=R2*X1
      S=X1
absg=ABS(x1)
!
!     INITIALIZE FOR COEFFICIENTS IN RECURSION FORMULAE
      P1=FL*E2
      P2=P1
      Q1=-E2
!
!     INITIALIZE ALPHA AND BETA
      ALP1=1.
      ALP2=1.+P2
      BET1=0.
      BET2=Q1
!
      IF(LL.EQ.0)GOTO 20
!
!     S AND SP FOR N=1
      X3=X3+2.
      X2=X2-1.
      X1=X1/X2
      SP=SP+X3*X1
      X1=R2*X1
      S=S+X1
absg=absg+ABS(x1)
!
!     LOOP FOR N=2 TO 2*LL
      DO 10 N=2,L2
!
!     CONTINUE CALCULATION OF FA AND PSI
      FN=N
      FA=FN*FA
      PS=PS+1./FN
!
!     CONTINUE CALCULATION OF S AND SP
      X3=X3+2.
      X2=X2-1.
      X1=X1/(X2*FN)
      SP=SP+X3*X1*ALP2
      X1=R2*X1
      S=S+X1*ALP2
absg=absg+ABS(x1*alp2)
!
!     COMPUTE COEFFICIENTS IN RECURSION FORMULAE
      P1=P1-E2
      P2=P2+P1
      Q1=Q1-E2
!     NOW HAVE P2=-N*(N-2*LL-1)*EPS/4
!     AND Q1=-N*EPS/2
!
!     NEW ALPHA AND BETA
      ALP0=ALP1
      ALP1=ALP2
      ALP2=ALP1+P2*ALP0
      BET0=BET1
      BET1=BET2
   10 BET2=BET1+P2*BET0+Q1*ALP0
!
!     NORMALIZE S AND SP, COMPLETE CALCULATION OF FA AND PS
      S=S*FA
absg=absg*ABS(fa)
      SP=SP*FA
      FA=FL2P1*FA
      PS=PS+1./FL2P1
!
!     COMPLETE CALCULATION OF ALPHA AND BETA
      P1=P1-E2
      P2=P2+P1
      Q1=Q1-E2
      ALP0=ALP1
      ALP1=ALP2
      BET0=BET1
      BET1=BET2
      BET2=BET1+P2*BET0+Q1*ALP0
!
   20 CONTINUE
!     NOW HAVE ALP1=ALPHA(2*LL+1)
!     AND BET1=BETA(2*LL+1), BET2=BETA(2*LL+2)
!
!     VALUE OF A=A(EPS,LL)
      A=ALP1
      A4=4.*A
      CL=2.*A*LOG(ABS(R2))
!D    FOR SINGLE PRECISION REPLACE DLOG BY ALOG AND DABS BY ABS
      CLP=2.*A/RHO
!
!  CALCULATE F AND FP AND CONTINUE CALCULATION OF S AND SP
!
!     CALCULATE A0,A1,D0,D1
      A0=(2.**LP1)/FA
      A1=-A0/FLP1
      PS=2.*PS*A
      D0=(BET1-PS)*A0
      D1=(BET2-PS-(2.+1./FLP1)*A)*A1
!
!     INITIALIZE F,FP, CONTINUE CALCULATION OF S,SP
!          -VALUES FOR N=0
!           U0 AND V0
      FNPLP1=FLP1
      C1=RHO**LL
      U0=A0*C1
      V0=D0*C1
      FP=FNPLP1*U0
      SP=SP+FNPLP1*V0
      U0=U0*RHO
      V0=V0*RHO
      F=U0
absf=ABS(f)
      S=S+V0
absg=absg+ABS(v0)
      W1=F*(CLP*F+SP)-FP*S
      NNN=0
!
!          - VALUES FOR N=1
!            U1 AND V1
      FNPLP1=FNPLP1+1.
      C1=C1*RHO
      U1=A1*C1
      V1=D1*C1
      FP=FP+FNPLP1*U1
      SP=SP+FNPLP1*V1
      U1=U1*RHO
      V1=V1*RHO
      F=F+U1
absf=absf+ABS(u1)
      S=S+V1
absg=absg+ABS(v1)
      W2=F*(CLP*F+SP)-FP*S
      DW2=ABS(W2-W1)
!
!     INITIALIZE FOR COEFFICIENTS IN RECURSION FORMULAE
      P1=-2.*FLP1
      P2=P1
      Q1=A4+2.*A*FL2P1
      REPS=RHO*EPS
!
!     LOOP FOR N=2 TO 100
      DO 40 N=2,100
!
!     COMPUTE COEFFICIENTS IN RECURSION FORMULAE
      P1=P1-2.
      P2=P2+P1
      Q1=Q1+DBLE(A4)
!     NOW HAVE P2=-N*(N+2*LL+1)
!     AND Q1=2*A*(2*N+2*LL+1)
!
!      COMPUTE U2=U(N) AND V2=V(N)
      U2=(2.*U1+DBLE(REPS)*U0)/DBLE(P2)
      V2=(2.*V1+DBLE(REPS)*V0+Q1*U2)/DBLE(P2)
!
!     INCREMENT FP AND SP
      FNPLP1=FNPLP1+1.
      FP=FP+DBLE(FNPLP1)*U2
      SP=SP+DBLE(FNPLP1)*V2
!
!     INCREMENT F AND S
      U2=U2*DBLE(RHO)
      V2=V2*DBLE(RHO)
      F=F+U2
absf=absf+ABS(u2)
      S=S+V2
absg=absg+ABS(v2)
!
!     CALCULATE WRONSKIAN
      W1=W2
      DW1=DW2
      W2=F*(CLP*F+SP)-FP*S
      DW2=ABS(W2-W1)
!
!     CONVERGENCE TEST
      K=N+1
      IF(DW1.GT.ACC2)GOTO 30
      IF(DW2.GT.ACC2)GOTO 30
      GOTO 50
!
!  NEW U0,U1,V0,V1
   30 U0=U1
      U1=U2
      V0=V1
      V1=V2
!
   40 CONTINUE
!
!  NOT CONVERGED
!
      IERR=2
      ACTACC=ABS(0.25*W2-1.)
      GOTO 60
!
!  CONVERGED
!
   50 ACTACC=ABS(0.25*W2-1.)
      IF(ACTACC.GT.ACC)IERR=1
!
!  COMPLETE CALCULATION OF G AND GP
!
   60 absg=absg*r2pi
      G=(S+CL*F)*R2PI
      GP=(SP+CL*FP+CLP*F)*R2PI

!
      RETURN
!
      END
!
!**********************************************************
!
      FUNCTION FKHI(E,L,AC)
use nrtype
!implicit real(kind=dbl)(a-h,o-z)
!
!  CALCULATES REAL PART OF (PSI(L+1+I*GAM)+PSI(GAM-L) - 2*LN(GAM))
!  WHERE E = 1/(GAM**2).
!  THIS IS REQUIRED FOR CALCULATION OF SCRIPT G.
!
implicit real(kind=dbl)(a-h,o-y)
      IMPLICIT COMPLEX(Z)
      PARAMETER(ONE=1.0d0,P0=ONE/252.d0)
!
      FKHI=0.
      IF(E.EQ.0)RETURN
!
      AC1=(20.*AC)**.333
!
      IF(E.GT.AC1)GOTO 100
!
      C=0.
      IF(L.EQ.0)GOTO 20
      A1=1.
      A2=-E
      A3=E+E
      DO 10 I=1,L
      A2=A2+A3
      A1=A1+A2
   10 C=C+FLOAT(I)/A1
   20 FKHI=E*((((1.05*E+1.)*E+2.1)*E+21.)*P0+C)
      RETURN
!
  100 AC1=1./SQRT(AC1)
      FKHI=0.
      FL=FLOAT(L+1)
      IF(FL.GT.AC1)GOTO 300
!
      N=AC1
      FL=N+1
      L1=L+1
      DO 210 I=L1,N
      FI=I
  210 FKHI=FKHI+FI/(1.+E*FI*FI)
      FKHI=-FKHI*E
!
  300 X1=FL*E
      X=1.+X1*FL
      ZE=CMPLX(FL,1./SQRT(E))
      ZE=-1./(ZE*ZE)
      FKHI=FKHI+.5*(LOG(X)-(X1/X))+REAL((((1.05*ZE+1.)*ZE&
     & +2.1)*ZE+21.)*ZE)*P0
!
      RETURN
      END
!***************************************************************
!
      SUBROUTINE ABG(E,L,AC,A,BG)
use nrtype
implicit real(kind=dbl)(a-h,o-z)
!
!  COMPUTES FUNCTION G(KAPPA,L) TO ACCURACY AC FOR E<O.
!
!      PARAMETER(PI=3.141592653589793238462643)
!
      X=1./SQRT(-E)
!
!  CALCULATION OF A AND EAC=E*A*C
      IF(L.GT.0)GOTO 2
      A=1.
      EAC=0.
      GOTO 20
    2 IF(L.GT.1)GOTO 4
      A=1.+E
      EAC=E
      GOTO 20
    4 IF(X.LT.FLOAT(L+1)) GOTO 12
      C=0.
      A=1.
      A1=1
      A2=-E
      A3=2.*E
      DO 10 I=1,L
      A2=A2+A3
      A1=A1+A2
      A=A*A1
   10 C=C+FLOAT(I)/A1
      EAC=E*A*C
      GOTO 20
!  CASE OF X.LT.(L+1)
   12 A=1.
      A1=0.
      DO 16  I=1,L
      A=A*(1.+FLOAT(I*I)*E)
      A2=FLOAT(I)
      DO 14 J=1,L
      IF(J.EQ.I) GOTO 14
      A2=A2*(1.+FLOAT(J*J)*E)
   14 CONTINUE
   16 A1=A1+A2
      EAC=E*A1
!
!  COMPUTE A1=PI*BG/A-E*C=1/(2*X)+PSI(X)-LN(X)
   20 A1=0.
!  TEST CONVERGENCE OF ASYMPTOTIC EXPANSION
      XN=(754.*AC)**(-.125)
      IF(X.GT.XN)GOTO 40
!  USE RECURRENCE FORMULAE
      N=XN-X+1
      XN=X+N
      Ee=-1./(XN*XN)
      A1=A1-.5*(1./X+1./XN)+LOG(XN/X)
      IF(N.LT.2)GOTO 40
      N=N-1
      DO 30 I=1,N
   30 A1=A1-1./(X+FLOAT(I))
!  USE ASYMPTOTIC EXPANSION
   40 A1=A1+(((1.05*Ee+1.)*Ee+2.1)*Ee+21)*Ee*.003968253968
!
!  COMPLETE CALCULATION
      BG=(A*A1+EAC)/PI
      RETURN
      END
!c************************************
!*****  input data fgh.in
! f
!0 10 0.5 1.e-10

!*****  output file fgh.out
! L= 0, E= 1.00000000E-02, R=-5.00000000E+00, AC= 1.00E-10
! f=-2.57088551E+02 [ 0.0], fp= 1.75566492E+02
! g= 6.63000089E-02 [ 3.5], gp=-4.77527285E-02
! h= 1.96337031E-03 [ 5.0], hp= 1.13547543E-03

! L= 0, ETA= 1.00000000E+01, RHO= 5.00000000E-01, AC= 1.00E-10
! F = 2.31408409E-12 [ 0.0], FP= 1.58029451E-11
! G = 3.42630245E+10 [ 5.0], GP=-1.98153259E+11
!                                                                            ****
