!!!*************************************************************
! 文件/File: pleg.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: pleg.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


function pleg(l,mm,theta)
  use nrtype,only : dbl,i4b,Pi
  implicit none
  integer(kind=i4b) :: l,m,mm
  real(kind=dbl) :: x,y,plgndr,sphe_arm_s,sphe_arm_c,pleg,norm,factrl,theta,fact1,fact2
  external plgndr,factrl

  if (mm.lt.0) then
     m=abs(mm)
  else
     m=mm
  end if
  x= cos(theta)

  fact1=factrl(l-m)
  fact2=factrl(l+m)
  norm=sqrt(dble((2*l+1))*fact1/(4.d0*Pi*fact2))
  pleg=plgndr(l,m,x)*norm

  if (mm.lt.0) then
     pleg=pleg*(-1)**m
  end if

end function pleg

function der_plg(l,m,theta)
  use nrtype,only : dbl,i4b,Pi
  use Brep, only : epsilon_custom
  implicit none
  real(kind=dbl) :: der_plg
  real(kind=dbl) :: x,plgndr,fact1,fact2,theta,factrl,norm
  real(kind=dbl), allocatable, save :: der_plg_lim(:,:,:)
  integer(kind=i4b) l,m,abs_m,mm,ll
  integer(kind=i4b),save :: counter=0,nx=11,ny=12,lmax=40,mmax
  external plgndr,factrl
  counter=counter+1
  abs_m=abs(m)
  if (counter==1) then
     mmax=lmax
     allocate(der_plg_lim(0:lmax,0:mmax,2))
     open(unit=nx,file='limit_values.1',status='old',action='read')
     open(unit=ny,file='limit_values.2',status='old',action='read')
     do ll=0,lmax
        do mm=0,ll
           read(ny,*)der_plg_lim(ll,mm,1)
           read(nx,*)der_plg_lim(ll,mm,2)
        end do
     end do
     close(nx)
     close(ny)
  end if
  write(6,*)'after read',theta,l,m
  x=cos(theta)
! Limit at 0, Pi of der_pleg is done reading in a table made in Mathematica
! analytically
! Need to factor out
! 2/sqrt(2),multiplied in parent subroutine
! Need to multiply limit at Pi by -1**m, it is wrong in Mathematica 
  if (theta.lt.0.d0+epsilon_custom) then
     der_plg=der_plg_lim(l,abs_m,1)*sqrt(2.d0)/2.d0 
  else if (theta.gt.Pi-epsilon_custom) then
     der_plg=der_plg_lim(l,abs_m,2)*sqrt(2.d0)/2.d0*(-1)**abs_m
  else
     if ((l-1).ge.m) then
        der_plg=(l*x*plgndr(l,m,x)-(l+m)*plgndr(l-1,m,x))/(sin(theta))
     else
        der_plg=(l*x*plgndr(l,m,x))/(sin(theta))
     end if
     fact1=factrl(l-m)
     fact2=factrl(l+m)
     norm=sqrt((2*l+1)*fact1/(4.d0*Pi*fact2))
     der_plg=der_plg*norm
  end if
  write(6,*)der_plg
end function der_plg

function new_der_plg(l,m,x)
  use nrtype,only : dbl,i4b
  implicit none
  real(kind=dbl) new_der_plg
  real(kind=dbl),intent(in) :: x
  integer(kind=i4b),intent(in) :: l,m
  real(kind=dbl) :: plgndr
  external :: plgndr
  write(120,*)l,m,x
  if ((l-1).ge.m) then
     new_der_plg=(l*m*plgndr(l,m,x)-(l+m)*plgndr(l-1,m,x))/(x*x-1)
  else
     new_der_plg=(l*m*plgndr(l,m,x))/(x*x-1)
  end if
end function new_der_plg


FUNCTION plgndr(l,m,x)
! SUBROUTINE TAKEN FROM NUMERICAL RECIPES
  use nrtype,only : dbl,i4b
  implicit none 
  INTEGER(kind=i4b) l,m 
  REAL(KIND=DBL) plgndr,x
  !Computes the associated Legendre polynomial P ml (x). Here m and l are integers satisfying 0 <= m <= l, while x lies in the range -1 <= x <= 1. 
  INTEGER(kind=i4b) i,ll 
  REAL(KIND=DBL) fact,pll,pmm,pmmp1,somx2 
  if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)pause 'bad arguments in plgndr' 
  pmm=1. 
  !Compute P mm . 
  if(m.gt.0) then
     somx2=sqrt((1.0d0-x)*(1.0d0+x)) 
     fact=1. 
     do 11 i=1,m
        pmm=-pmm*fact*somx2 
        fact=fact+2.0d0 
11   enddo
  endif
  if(l.eq.m) then
     plgndr=pmm 
  else
     pmmp1=x*(2*m+1)*pmm 
     !Compute P mm+1. 
     if(l.eq.m+1) then
        plgndr=pmmp1 
     else 
        !Compute P ml , l > m + 1.
        do 12 ll=m+2,l
           pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m) 
           pmm=pmmp1 
           pmmp1=pll 
12      enddo
        plgndr=pll 
     endif
  endif
  !write(6,*)plgndr
  return 
END FUNCTION plgndr






FUNCTION factrl(n)
  use nrtype,only : dbl,i4b
  implicit none 
  INTEGER(kind=i4b) n 
  REAL(KIND=DBL) factrl 
  !USES gammln

  !Returns the value n! as a floating-point number. 
  INTEGER(kind=i4b) j,ntop 
  REAL(KIND=DBL) a(33),gammln 
  !Table to be filled in only as required. 
  SAVE ntop,a 
  DATA ntop,a(1)/0,1./ 
  !Table initialized with 0! only. 
  if (n.lt.0) then

     pause 'negative factorial in factrl' 
  else if (n.le.ntop) then 
     !Already in table.
     factrl=a(n+1) 
  else if (n.le.32) then 
     ! Fill in table up to desired value.
     do 11 j=ntop+1,n

        a(j+1)=j*a(j) 
11   enddo
     ntop=n 
     factrl=a(n+1) 
  else 
     !Larger value than size f table is required.Actually,this big a value is going to over  w n many computers,but no harm in trying. 
     !factrl=exp(gammln(n+1.))
  endif
  return 
END FUNCTION factrl


function f_phi(angle_phi,m,f)
  use nrtype,only : dbl,i4b
  implicit none
  real(kind=dbl) angle_phi,f_phi
  integer(kind=i4b) m,f
  !write(328,*)angle_phi,m,f
  if (f.eq.1) then
     f_phi=cos(m*angle_phi)
  else if (f.eq.2) then
     f_phi=sin(m*angle_phi)
  end if
  !write(331,*)angle_phi
end function f_phi



function d_phi(angle_phi,m,f)
  use nrtype,only : dbl,i4b
  implicit none
  real(kind=dbl) angle_phi,d_phi
  integer(kind=i4b) m,f
  if (f.eq.1) then
     d_phi=-m*sin(m*angle_phi)
  else if (f.eq.2) then
     d_phi=m*cos(angle_phi)
  end if

end function d_phi

FUNCTION gammaln(xx) 
  use nrtype,only : dbl,i4b

  real(kind=dbl) :: gammaln,xx 
  !Returns the value ln [ (xx ]for xx > 0 .
  INTEGER(kind=i4b) :: j 
  real(kind=dbl) :: ser,stp,tmp,x,y,cof(6) 
  !Internal arithmetic will be done in double precision,a nicety that you can mit if  ve- gure accuracy is god enough. 
  SAVE cof,stp 
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, -.5395239384953d-5,2.5066282746310005d0/ 
  x=xx 
  y=x 
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp 
  ser=1.000000000190015d0 
  do  j=1,6 
     y=y+1.d0
     ser=ser+cof(j)/y 
  enddo
  gammaln=tmp+log(stp*ser/x) 
  return 
END FUNCTION gammaln


!     ALGORTHM 404 COLLECTED ALGORITHMS FROM ACM.
!     ALGORITHM APPEARED IN COMM. ACM, VOL. 14, NO. 01,
!     P. 048.
FUNCTION CGAMMA(Z)
  use nrtype,only : dbl,i4b,dpc

  implicit none
  COMPLEX(kind=dpc):: Z,ZM,T,TT,SUM,TERM,DEN,CGAMMA,PI,A
  real(kind=dbl),DIMENSION(12):: C
  LOGICAL REFLEK
  real(kind=dbl) x,y,tol,xdist
  integer(kind=i4b) :: m,i,j,iout    
  ! SET IOUT FOR PROPER OUTPUT CHANNEL OF COMPUTER SYSTEM FOR
  ! ERROR MESSAGES
  IOUT = 3
  PI = (3.141592653589793238462643383279502884197d0,0.0d0)
  X = REAL(Z)
  Y = AIMAG(Z)
  ! TOL = LIMIT OF PRECISION OF COMPUTER SYSTEM IN SINGLE PRECISI
  TOL = 1.0d-7
  REFLEK = .TRUE.
  ! DETERMINE WHETHER Z IS TOO CLOSE TO A POLE
  ! CHECK WHETHER TOO CLOSE TO ORIGIN
  IF(X.GE.TOL) GO TO 20
  ! FIND THE NEAREST POLE AND COMPUTE DISTANCE TO IT
  XDIST = X-INT(X-.5)
  ZM = CMPLX(XDIST,Y)
  IF(ZABS(ZM).GE.TOL) GO TO 10
  ! IF Z IS TOO CLOSE TO A POLE, PRINT ERROR MESSAGE AND RETURN
  ! WITH CGAMMA = (1.E7,0.0E0)
  WRITE(IOUT,900) Z
  CGAMMA = (1.d7,0.d0)
  RETURN
  ! FOR REAL(Z) NEGATIVE EMPLOY THE REFLECTION FORMULA
  ! GAMMA(Z) = PI/(SIN(PI*Z)*GAMMA(1-Z))
  ! AND COMPUTE GAMMA(1-Z).  NOTE REFLEK IS A TAG TO INDICATE THA
  ! THIS RELATION MUST BE USED LATER.
10 IF(X.GE.0.0) GO TO 20
  REFLEK = .FALSE.
  Z = (1.0d0,0.0d0)-Z
  X = 1.0-X
  Y = -Y
  ! IF Z IS NOT TOO CLOSE TO A POLE, MAKE REAL(Z)>10 AND ARG(Z)<P
20 M = 0
40 IF(X.GE.10.) GO TO 50
  X = X + 1.0d0
  M = M + 1
  GO TO 40
50 IF(ABS(Y).LT.X) GO TO 60
  X = X + 1.0d0
  M = M + 1
  GO TO 50
60 T = CMPLX(X,Y)
  TT = T*T
  DEN = T
  ! COEFFICIENTS IN STIRLING*S APPROXIMATION FOR LN(GAMMA(T))
  C(1) = 1./12.d0
  C(2) = -1./360.d0
  C(3) = 1./1260.d0
  C(4) = -1./1680.d0
  C(5) = 1./1188.d0
  C(6) = -691./360360.d0
  C(7) = 1./156.d0
  C(8) = -3617./122400.d0
  C(9) = 43867./244188.d0
  C(10) = -174611./125400.d0
  C(11) = 77683./5796.d0
  SUM = (T-(.5d0,0.0d0))*ZLOG(T)-T+CMPLX(.5d0*LOG(2.d0*3.14159d0),0.0d0)
  J = 1
70 TERM = C(J)/DEN
  ! TEST REAL AND IMAGINARY PARTS OF LN(GAMMA(Z)) SEPARATELY FOR
  ! CONVERGENCE.  IF Z IS REAL SKIP IMAGINARY PART OF CHECK.
  IF(ABS(REAL(TERM)/REAL(SUM)).GE.TOL) GO TO 80
  IF(Y.EQ.0.0) GO TO 100
  IF(ABS(AIMAG(TERM)/AIMAG(SUM)).LT.TOL) GO TO 100
80 SUM = SUM + TERM
  J = J + 1
  DEN = DEN*TT
  ! TEST FOR NONCONVERGENCE
  IF(J.EQ.12) GO TO 90
  GO TO 70
  ! STIRLING*S SERIES DID NOT CONVERGE.  PRINT ERROR MESSAGE AND
  ! PROCEDE.
90 WRITE(IOUT,910) Z
  ! RECURSION RELATION USED TO OBTAIN LN(GAMMA(Z))
  ! LN(GAMMA(Z)) = LN(GAMMA(Z+M)/(Z*(Z+1)*...*(Z+M-1)))
  ! = LN(GAMMA(Z+M)-LN(Z)-LN(Z+1)-...-LN(Z+M
100 IF(M.EQ.0) GO TO 120
  DO I = 1,M
     A = CMPLX(I*1.-1.,0.0d0)
110  SUM = SUM-ZLOG(Z+A)
  end do
  ! CHECK TO SEE IF REFLECTION FORMULA SHOULD BE USED
120 IF(REFLEK) GO TO 130
  SUM = ZLOG(PI/ZSIN(PI*Z))-SUM
  Z = (1.0d0,0.0d0) -Z
130 CGAMMA = ZEXP(SUM)
  RETURN
900 FORMAT(1X,2E14.7,10X)
910 FORMAT(14X,4HZ = ,2E14.7)
END FUNCTION CGAMMA
