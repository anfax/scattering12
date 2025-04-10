!!!*************************************************************
! 文件/File: potential.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: potential.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


real(kind=8) function hara(density_cuberoot,energy,abs_E_last_bound)
  use nrtype, only : dbl,Pi
  use Calc_func, only : Energy_Hara
  implicit none
  real(kind=dbl) :: density_cuberoot,energy,abs_E_last_bound
  real(kind=dbl), parameter :: third=.33333333333333333333333d0,three_Pi_square_cuberoot=3.093667726280135930969d0,two_over_Pi=.63661977236758134307d0
  real(kind=dbl) :: eta,k_Fermi,kappa,k_hara,func_eta
  !abs_E_last_bound=abs(E_last_bound)
  !k_Fermi=(3.d0*Pi*Pi)**third*density_cuberoot
  !k_Fermi=(3.d0*Pi*Pi*density_cuberoot)**third
  k_Fermi=three_Pi_square_cuberoot*(density_cuberoot)
  !kappa=sqrt(2.d0*energy)
  !k_hara=sqrt(k_Fermi*k_Fermi+2.d0*abs_E_last_bound+kappa*kappa)
  !k_hara=sqrt(k_Fermi*k_Fermi+2.d0*abs_E_last_bound+2*energy)
  k_hara=sqrt(k_Fermi*k_Fermi+Energy_Hara)
  eta=k_hara/k_Fermi
  func_eta=.5d0+((1.0d0 -eta*eta)/(4.d0*eta))*log(abs((1.d0+eta)/(1.d0 -eta)))
  !write(308,*)func_eta
  !hara= - 2.d0 /Pi*k_Fermi*func_eta
  hara= - two_over_Pi*k_Fermi*func_eta
end function hara


!real(kind=dbl) function hfege(chg,xksq,xkoop)
!  use nrtype
!	Subroutine that calculates the Hara exchange potential (from SCATPAK)
!  implicit none
!!  real(kind=dbl) :: ve,v,chg,xksq,xkoop,hkmx,hkmxsq,hkcsq,hkc,arg,vhara
!  real(kind=dbl) :: third,piov4
!  data third/.33333333333333333d0/
!  piov4=Pi/4
!  hkmx=(3.*piov4*chg)**third
!  hkmxsq=hkmx*hkmx
!  hkcsq=xksq+xkoop+hkmxsq
!  hkc=dsqrt(hkcsq)
!  arg=(hkc-hkmx)/(hkc+hkmx)
!  arg=dlog(arg)
!     vhara=2.*((xksq+xkoop)*arg/(2.*hkc)+hkmx)/pi
!     ve=v-vhara
!!  hfege=2.*((xksq+xkoop)*arg/(2.*hkc)+hkmx)/pi  
!  return
!end function hfege


!subroutine charge(charge_density)
!	Calculates the corrections to the potential (Haara & Polarization)
!!  use nrtype
!  use Brep
!  use potential_data
!  use open_information, only : kappa
!  integer element
!  integer,parameter :: nx=13 
!  real(kind=dbl),dimension(n):: charge_density
!	unit name: nx=13 for the charge

!  call V_conversion(charge_density,n,nx)
!  kappasq=kappa*kappa
!  do element=1,n
!    U(element)=U(element)-hfege(charge_density(element),kappasq,hi_bound_state)
!    U(element)=U(element)+v_pol(alpha0,r_cutoff,r_average(element))
!  end do

!end subroutine charge



!real(kind=dbl) function v_pol(alpha0,r_cutoff,r)
!  use nrtype
!  real(kind=dbl) :: alpha0,r_cutoff,vp,v
!  r2=r*r
!  r4=r2*r2
!  r6=r4*r2
!  r_cutoff2=r_cutoff*r_cutoff
!  r_cutoff6=r_cutoff2*r_cutoff2*r_cutoff2
!  v_pol=alpha0/(r4)*(1-nep_e**(r6/r_cutoff6))
!end function v_pol
