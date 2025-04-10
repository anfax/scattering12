!!!*************************************************************
! 文件/File: Hydrogen_perturbed.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: Hydrogen_perturbed.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


module pass_quantum_numbers_H
  use nrtype
  integer(kind=dbl) :: n,l,m,np,lp,mp,NumSteps
  real(kind=dbl) :: Z,R0,DeltaE,Emin,cos_alpha,cos_beta,cos_gamma
end module pass_quantum_numbers_H

!program Hydrogen_perturbed
module Hydrogen
contains
subroutine Hydrogen_perturbed(dipole,H,eigen,dime)
  use nrtype,only : dbl,dpc,i4b,Pi
  use gauss2dint
  use pass_quantum_numbers_H, only:cos_alpha,cos_beta,cos_gamma
  implicit none
  integer(kind =i4b) :: n,l,m,nmax,dime,i,j,nn,counter
  real(kind=dbl) :: eigen(dime),H(dime,dime)
  real(kind=dbl),parameter :: Debye_to_au=2.541765d0
  real(kind=dbl) :: dipole(5),mu
  complex(kind=dpc), allocatable :: eigen_complex(:),H_complex(:,:) 

  ! Parameters
  if (dipole(1).ne.0) then
  dipole(1:4)=dipole(1:4)/Debye_to_au
  cos_alpha=dipole(2)/dipole(1)
  if (dipole(5).gt.0) then
    cos_beta=dipole(3)/dipole(1)
    cos_gamma=dipole(4)/dipole(1)
  else
    cos_beta=dipole(4)/dipole(1)
    cos_gamma=dipole(3)/dipole(1)
  end if
  else
   cos_alpha=0.d0;cos_beta=0.d0;cos_gamma=1.d0 !Case of zero dipole
  end if
  !call read_input(nmax)
  write(6,*)'dipole moment=',dipole
  write(6,*)'direction cosines=',cos_alpha,cos_beta,cos_gamma

  !  nmax=3
  !dime=0
  ! parameters
  !dime=(nmax+1)**2
  nmax=anint(sqrt(dble(dime))-1)
  write(6,*)'dimension=',dime,nmax
  !allocate(eigen(dime))
  !allocate(H(dime,dime))

  call CalcHamiltonian(nmax,dime,H,dipole)

  !write(36,*)H(:,:)
  !allocate(H_complex(dime,dime))
  allocate(eigen_complex(dime))
  !H_complex=cmplx(H,0.d0)
  !eigen_complex=cmplx(0.d0,0.d0)
  !call diag_complex(H_complex,dime,dime,eigen_complex)
  !write(23,*)eigen_complex
  call diag_symm(H,dime,dime,eigen)

  eigen_complex(:)=(-1.d0+sqrt(1.d0+4.d0*cmplx(eigen(:),0.d0)))/2.d0
  !--------------------------
  ! Outputs
  write(24,*)eigen_complex(:)
  do i=1,dime
  write(22,*)eigen(i),dipole(1)
  end do
  write(25,*)H(:,1)
  !--------------------------
  !call continuum_main(H,dime,dime,eigen)
  !deallocate(H_complex)
  deallocate(eigen_complex)
  !deallocate(H)
  !deallocate(eigen)
end subroutine Hydrogen_perturbed

subroutine CalcHamiltonian(nmax,dime,H,dipole)
  use nrtype,only : dbl,dpc,i4b,Pi
  !use gauss3dint
  use pass_quantum_numbers_H
  use gauss2dint 
  implicit none
  integer(kind=i4b),intent(in) :: nmax,dime
  integer(kind =i4b) :: index,indexx,num_pts,lmax,ngp
  real(kind=dbl) :: ss,func_integrate,a,b,r_lo,r_hi,theta_lo,theta_hi,phi_lo,phi_hi,func_theta_H,&
       & func_r,s_r,s_theta,s_theta_diag,func_phi_H,s_phi
  real(kind=dbl),allocatable :: weigh(:),xabsc(:)
  real(kind=dbl),intent(in) :: dipole(5)
  real(kind=dbl),intent(out) :: H(dime,dime)
  external :: func_theta_H,func_theta_diag_H,func_phi_H,func_hydrogen_2D
  num_pts=30
  ngp=num_pts
  r_lo=0.d0
  r_hi=200.d0
  theta_lo=0.d0
  theta_hi=Pi
  phi_lo=0.d0
  phi_hi=2.d0*Pi
  index=0
  ss=0.d0
  s_r=0.d0
  s_theta=0.d0
  s_phi=0.d0
  lmax=nmax
  allocate(weigh(ngp))
  allocate(xabsc(ngp))

  call gauleg(ngp,xabsc,weigh)
  write(6,*)'lmax=',lmax
  do l=0,lmax
     do m=-l,l
        index=index+1
        indexx=0
        do lp=0,lmax
           do mp=-lp,lp
              indexx=indexx+1
              if (index.eq.indexx) then
!                 call qgaus(func_theta_H,theta_lo,theta_hi,s_theta)
                 call qgaus(func_theta_diag_H,theta_lo,theta_hi,s_theta_diag)
!,xabsc,weigh,ngp)
!                 call qgaus(func_phi_H,phi_lo,phi_hi,s_phi)
!,xabsc,weigh,ngp)
!----------------------------------------
!$ 2nd term is always zero, eliminated
                 H(index,indexx)=(l*(l+1))*s_theta_diag-0.d0*(2.d0*dipole(1))*s_theta*s_phi
!-----------------------------------------
                 !write(12,*)index,n,l,m,np,lp,mp,s_theta,s_theta_diag
              else
                 !if (m.eq.mp) then
                    !                       call qgaus(func_r,r_lo,r_hi,s_r)
!                    call qgaus(func_theta_H,theta_lo,theta_hi,s_theta)
!,xabsc,weigh,ngp)
!                    call qgaus(func_phi_H,phi_lo,phi_hi,s_phi)
!,xabsc,weigh,ngp)
                    counter=0
                    ss=qgss2d(func_hydrogen_2D,0.d0,Pi,0.d0,2.d0*Pi,ngp)
                    H(index,indexx)=-2.d0*ss*dipole(1)
!s_theta*s_phi
                    write(37,*)ss*dipole(1),l,lp
                    !write(13,*)n,l,m,np,lp,mp,s_r,ss
                    !write(14,*)n-l-1,2*l+1,m,np-lp-1,2*lp+1,mp,s_theta
                    !s_r*s_theta*dipole
                 !else
                 !   H(index,indexx)=0.d0
                 !end if
              end if
           end do
        end do
     end do
  end do
  deallocate(weigh)
  deallocate(xabsc)

  write(6,*)'CalcHamiltonian end'
end subroutine CalcHamiltonian

function func_integrate_H(rr,theta,phi)
  use nrtype, only : dbl,i4b
  use gauss3dint
  use pass_quantum_numbers_H
  implicit none
  real(kind=dbl) :: func_integrate_H,lag_pol(0:n),rescaled_variable,factrl,pleg&
       & ,norm,normp
  real(kind=dbl),intent(in) :: rr,theta,phi
  external :: factrl,pleg
  rescaled_variable=2.d0*rr/n
  norm=Sqrt((2.d0/n)**3*factrl(n-l-1)/(2*n*(factrl(n+l))**3))*exp(-rr/n)*(2.d0*rr/n)**l
  normp=Sqrt((2.d0/np)**3*factrl(np-lp-1)/(2*np*(factrl(np+lp))**3))*exp(-rr/np)*(2.d0*rr/np)**lp
  !call laguerre_lnm(n-l-1,2*l+1,rescaled_variable,lag_pol)
  !call laguerre_lnm(np-lp-1,2*lp+1,rescaled_variable,lag_pol)
  func_integrate_H=pleg(l,m,theta)*pleg(lp,mp,theta)*sin(theta)
end function func_integrate_H




function func_r(rr)
  use nrtype, only : dbl,i4b,Pi
  use pass_quantum_numbers_H
  implicit none
  real(kind=dbl) :: func_r,lag_pol(0:n+l),lag_polp(0:np+lp),factrl,rescaled_variable,norm,normp,&
       & rescaled_variablep
  real(kind=dbl),intent(in) :: rr
  external :: factrl,pleg
  rescaled_variable=2.d0*rr/n
  rescaled_variablep=2.d0*rr/np
  norm=Sqrt((2.d0/n)**3 * factrl(n-l-1)/(2*n*(factrl(n+l)))) * exp(-rr/n)*(2.d0*rr/n)**l 
  normp=Sqrt((2.d0/np)**3 * factrl(np-lp-1)/(2*np*(factrl(np+lp)))) * exp(-rr/np)*(2.d0*rr/np)**lp 
  !call laguerre_lnm(n-l-1,2*l+1,rescaled_variable,lag_pol)
  !call laguerre_lnm(np-lp-1,2*lp+1,rescaled_variablep,lag_polp)
  func_r=lag_pol(n-l-1)*lag_polp(np-lp-1)*norm*normp*(1.d0-exp(-(rr/2.d0)**6))
  !*rr**2

  !norm=Sqrt((2.d0/n)**3 * factrl(n-l-1)/(2*n*(factrl(n+l))**3)) *exp(-rescaled_variable/2.d0)*rescaled_variable**l
  !normp=Sqrt((2.d0/np)**3 * factrl(np-lp-1)/(2*np*(factrl(np+lp))**3)) *exp(-rescaled_variablep/2.d0)*rescaled_variablep**lp
  !call laguerre_lnm(n+l,2*l+1,rescaled_variable,lag_pol)
  !call laguerre_lnm(np+lp,2*lp+1,rescaled_variablep,lag_pol)
  !write(3002,*)lag_pol(n-l-1),n,rescaled_variable,n-l-1,2*l+1
  !func_r=-lag_pol(n+l)*lag_pol(np+lp)*norm*normp*rr**2
end function func_r



subroutine diag_symm(a,n,lda,eigen)
  use nrtype, only : dbl,i4b
  implicit none
  integer(kind=i4b) itype,lda,ldb,info,lwork,i,j,ldvr,ldvl,nrhs,n,m,counter2
  character jobvl,jobvr,jobz,uplo
  real(kind=dbl),dimension(n,n) :: A
  real(kind=dbl),dimension(:,:),allocatable ::vl,vr
  real(kind=dbl),dimension(:),allocatable :: work,alphar,alphai,den,wr,wi,w
  real(kind=dbl),dimension(n) :: eigen
  counter2=counter2+1
  uplo='u'
  jobz='v'
  jobvl='n'
  jobvr='v'
  lda=n
  ldb=lda
  ldvr=n
  ldvl=n
  lwork=3*n-1
  itype=1
  m=lda
  ldb=lda
  allocate(w(n))
  allocate(work(lwork))
  write(6,*)'before dsygv',jobz,uplo,n,lda
  call dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
  write(6,*)'sfter dsygv info=',info

  eigen=w
90 format(4d14.7)
  write(23,*)w
  write(6,*)'info=',info
  deallocate(work)
  deallocate(w)

end subroutine diag_symm

subroutine diag_complex(a,n,lda,eigen)
  use nrtype, only : dbl,i4b
  implicit none
  integer(kind=i4b) itype,lda,ldb,info,lwork,i,j,ldvr,ldvl,nrhs,n,m,counter2
  character jobvl,jobvr,jobz,uplo
  complex*16,dimension(n,n) :: A
  complex*16,dimension(:,:),allocatable ::vl,vr
  complex*16,dimension(:),allocatable :: work,wr,wi,w
  real(kind=dbl),dimension(:),allocatable :: alphar,alphai,den,rwork
  complex*16,dimension(n) :: eigen
  counter2=counter2+1
  uplo='u'
  jobz='v'
  jobvl='n'
  jobvr='v'
  lda=max(1,n)
  ldb=lda
  ldvr=n
  ldvl=n
  lwork=max(1,2*n)
  itype=1
  m=lda
  ldb=lda
  allocate(w(n))
  allocate(work(lwork))
  allocate(vr(n,n))
  allocate(vl(n,n))
  allocate(alphai(n))
  allocate(alphar(n))
  allocate(den(n))
  allocate(wr(n))
  allocate(wi(n))
  allocate(rwork(lwork))
  call ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK,LWORK, RWORK, INFO )
  if (info.ne.0) then
     write(6,*)'info=',info
  end if
  eigen=w
  deallocate(work)
  deallocate(vr)
  deallocate(vl)
  deallocate(alphai)
  deallocate(alphar)
  deallocate(den)
end subroutine diag_complex



SUBROUTINE qgaus(func,a,b,ss) 
use nrtype
REAL(kind=dbl) a,b,ss,func 
EXTERNAL func 
!Returns as ss the integral of the function func between a and b, by ten-point Gauss- Legendre
!integration: the function is evaluated exactly ten times at interior points in
!the range of integration. 
REAL(kind=dbl) dx,xm,xr,w(40),x(40)
INTEGER j
 SAVE w,x 

!data x/  0.38772417506050821933193444024624E-01&
!&,  0.11608407067525520848345128440802E+00&
!&,  0.19269758070137109971551685206515E+00&
!&,  0.26815218500725368114118434480860E+00&
!&,  0.34199409082575847300749248117919E+00&
!&,  0.41377920437160500152487974580371E+00&
!&,  0.48307580168617871290856657424482E+00&
!&,  0.54946712509512820207593130552952E+00&
!&,  0.61255388966798023795261245023069E+00&
!&,  0.67195668461417954837935451496149E+00&
!&,  0.72731825518992710328099645175493E+00&
!&,  0.77830565142651938769497154550650E+00&
!&,  0.82461223083331166319632023066610E+00&
!&,  0.86595950321225950382078180835462E+00&
!&,  0.90209880696887429672825333086849E+00&
!&,  0.93281280827867653336085216684521E+00&
!&,  0.95791681921379165580454099945276E+00&
!&,  0.97725994998377426266337028371290E+00&
!&,  0.99072623869945700645305435222137E+00&
!&,  0.99823770971055920034962270242059E+00/

!data w/0.77505947978424811263723962958326E-01&
!&,0.77039818164247965588307534283811E-01&
!&,0.76110361900626242371558075922495E-01&
!&,0.74723169057968264200189336261325E-01&
!&,0.72886582395804059060510683442519E-01&
!&,0.70611647391286779695483630855287E-01&
!&,0.67912045815233903825690108231923E-01&
!&,0.64804013456601038074554529566753E-01&
!&,0.61306242492928939166537996408400E-01&
!&,0.57439769099391551366617730910427E-01&
!&,0.53227846983936824354996479772260E-01&
!&,0.48695807635072232061434160448145E-01&
!&,0.43870908185673271991674686041717E-01&
!&,0.38782167974472017639972031290446E-01&
!&,0.33460195282547847392678183086411E-01&
!&,0.27937006980023401098489157507721E-01&
!&,0.22245849194166957261504324184209E-01&
!&,0.16421058381907888712863484882365E-01&
!&,0.10498284531152813614742171067278E-01&
!&,0.45212770985331912584717328781864E-02/

data x/&
&0.1951138325679399765435123D-01&
&,0.5850443715242066862899332D-01&
&,0.9740839844158459906327845D-01&
&,0.1361640228091438865592411D+00&
&,0.1747122918326468125593390D+00&
&,0.2129945028576661325723885D+00&
&,0.2509523583922721204931588D+00&
&,0.2885280548845118531091393D+00&
&,0.3256643707477019146191129D+00&
&,0.3623047534994873156190433D+00&
&,0.3983934058819692270243796D+00&
&,0.4338753708317560930623867D+00&
&,0.4686966151705444770360784D+00&
&,0.5028041118887849875936728D+00&
&,0.5361459208971319320198573D+00&
&,0.5686712681227097847254858D+00&
&,0.6003306228297517431547463D+00&
&,0.6310757730468719662479284D+00&
&,0.6608598989861198017359671D+00&
&,0.6896376443420276007712076D+00&
&,0.7173651853620998802540683D+00&
&,0.7440002975835972723165405D+00&
&,0.7695024201350413738656161D+00&
&,0.7938327175046054499486393D+00&
&,0.8169541386814634703711250D+00&
&,0.8388314735802552756166230D+00&
&,0.8594314066631110969771921D+00&
&,0.8787225676782138287037733D+00&
&,0.8966755794387706831943241D+00&
&,0.9132631025717576541647337D+00&
&,0.9284598771724457959530460D+00&
&,0.9422427613098726747522660D+00&
&,0.9545907663436349054934815D+00&
&,0.9654850890437992514522732D+00&
&,0.9749091405857277933856452D+00&
&,0.9828485727386290704182880D+00&
&,0.9892913024997555310265032D+00&
&,0.9942275409656882778920635D+00&
&,0.9976498643982376888994942D+00&
&,0.9995538226516306298800805D+00/













data w/&
&0.3901781365630665481128044D-01&
&,0.3895839596276953119862552D-01&
&,0.3883965105905196893177418D-01&
&,0.3866175977407646332707711D-01&
&,0.3842499300695942318521244D-01&
&,0.3812971131447763834420679D-01&
&,0.3777636436200139748977498D-01&
&,0.3736549023873049002670538D-01&
&,0.3689771463827600883915100D-01&
&,0.3637374990583597804396499D-01&
&,0.3579439395341605460286159D-01&
&,0.3516052904474759349552659D-01&
&,0.3447312045175392879436423D-01&
&,0.3373321498461152281667516D-01&
&,0.3294193939764540138283618D-01&
&,0.3210049867348777314805649D-01&
&,0.3121017418811470164244287D-01&
&,0.3027232175955798066122001D-01&
&,0.2928836958326784769276759D-01&
&,0.2825981605727686239675320D-01&
&,0.2718822750048638067441871D-01&
&,0.2607523576756511790296874D-01&
&,0.2492253576411549110511785D-01&
&,0.2373188286593010129319252D-01&
&,0.2250509024633246192622159D-01&
&,0.2124402611578200638871074D-01&
&,0.1995061087814199892889193D-01&
&,0.1862681420829903142873541D-01&
&,0.1727465205626930635858421D-01&
&,0.1589618358372568804490291D-01&
&,0.1449350804050907611696207D-01&
&,0.1306876159240133929378683D-01&
&,0.1162411412079782691646677D-01&
&,0.1016176604110306452083185D-01&
&,0.8683945269260858426409452D-02&
&,0.7192904768117312752675571D-02&
&,0.5690922451403198649269107D-02&
&,0.4180313124694895236739304D-02&
&,0.2663533589512681669293535D-02&
&,0.1144950003186941534544170D-02/









!  data w/.295524224714753,.269266719309996,.219086362515982,.149451349150581,.066671344308688/
!  data x/.148874338981631,.433339534129247,.679409568299024,.865063366688985,.973906528517172/
xm=0.5*(b+a)
xr=0.5*(b-a) 
ss=0
do  j=1,40
dx=xr*x(j)
ss=ss+w(j)*(func(xm+dx)+func(xm-dx)) 
enddo 
 ss=xr*ss 
return 
END subroutine qgaus
end module Hydrogen

function func_theta_H(theta)
  use nrtype, only : dbl,i4b,Pi
  use pass_quantum_numbers_H
  implicit none
  real(kind=dbl) :: func_theta_H,pleg,factrl
  real(kind=dbl),intent(in) :: theta
  external :: factrl,pleg
  func_theta_H=pleg(l,m,theta)*pleg(lp,mp,theta)*sin(theta)*sin(theta)

end function func_theta_H

function func_theta_diag_H(theta)
  use nrtype, only : dbl,i4b,Pi
  use pass_quantum_numbers_H
  implicit none
  real(kind=dbl) :: func_theta_diag_H,pleg,factrl
  real(kind=dbl),intent(in) :: theta
  external :: factrl,pleg
  func_theta_diag_H=pleg(l,m,theta)*pleg(lp,mp,theta)*sin(theta)*2.d0*Pi
end function func_theta_diag_H

function func_phi_H(phi)
  use nrtype, only : dbl,i4b,Pi
  use pass_quantum_numbers_H
  implicit none
  real(kind=dbl) :: func_phi_H,pleg,factrl,ff,ffp
  real(kind=dbl),intent(in) :: phi
  external :: factrl,pleg
!if (mod(abs(m),2).eq.0) then
  if (m.gt.0) then
    ff=cos((phi)*(abs(m)))*sqrt(2.d0)
  else if (m.lt.0) then
    ff=sin((phi)*(abs(m)))*sqrt(2.d0)
  else
    ff=1.d0
  end if
!else
!  if (m.gt.0) then
!    ff=sin((phi)*(abs(m)))*sqrt(2.d0)
!  else if (m.lt.0) then
!    ff=cos((phi)*(abs(m)))*sqrt(2.d0)
!  else
!    ff=1.d0
!  end if
!end if

!if (mod(abs(mp),2).eq.0) then
  if (mp.gt.0) then
    ffp=cos((phi)*(abs(mp)))*sqrt(2.d0)
  else if (mp.lt.0) then
    ffp=sin((phi)*(abs(mp)))*sqrt(2.d0)
  else
    ffp=1.d0
  end if
!else
!  if (mp.gt.0) then
!    ffp=sin((phi)*(abs(mp)))*sqrt(2.d0)
!  else if (mp.lt.0) then
!    ffp=cos((phi)*(abs(mp)))*sqrt(2.d0)
!  else
!    ffp=1.d0
!  end if
!end if
  func_phi_H=ff*ffp*Sqrt(2.d0)/2.d0*(-cos(phi+0.d0*Pi/4.d0)+sin(phi))
end function func_phi_H



function func_hydrogen_2D(theta,phi)
  use nrtype, only : dbl,i4b,Pi
  use pass_quantum_numbers_H
  implicit none
  real(kind=dbl) :: func_hydrogen_2D,pleg,factrl,ff,ffp
  real(kind=dbl),intent(in) :: theta,phi
  external :: factrl,pleg
  if (m.gt.0) then
    ff=cos((phi)*(abs(m)))*sqrt(2.d0)
  else if (m.lt.0) then
    ff=sin((phi)*(abs(m)))*sqrt(2.d0)
  else
    ff=1.d0
  end if
  if (mp.gt.0) then
    ffp=cos((phi)*(abs(mp)))*sqrt(2.d0)
  else if (mp.lt.0) then
    ffp=sin((phi)*(abs(mp)))*sqrt(2.d0)
  else
    ffp=1.d0
  end if
  func_hydrogen_2D=ff*ffp*pleg(l,m,theta)*pleg(lp,mp,theta)*sin(theta)*&
       & (cos_alpha*sin(theta)*cos(phi)+cos_beta*sin(theta)*sin(phi)+cos_gamma*cos(theta))*(-1)**abs(mp)
end function func_hydrogen_2D

