!!!*************************************************************
! 文件/File: separate_integral.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: separate_integral.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


function separate_integral(a,b,num_pts,element_local)
  use nrtype, only : Pi,i4b,dp
  use Calc_func, only : an_r,an_theta,an_phi,der_r,bas_r,k_ri,k_rj,k_thetai,k_thetaj,k_phii,k_phij
  !use control
  !use open_information
  use function_integr_tab, only : der_tab,bas_tab,angle_theta_tab,radius_tab
  implicit none
  integer(kind=i4b) :: ii, count_localx,count_localy,count_localz
  integer(kind=i4b) :: num_pts,element_local
  real(kind=dp), intent(in) :: a,b
  real(kind=dp) :: x,y,zeta,func_der_r,func_der_theta,func_der_phi,func_1_r,&
       & func_1_theta,func_1_phi,func_2_theta,separate_integral,&
       & rad,sin_theta,pm,pl
  real(kind=dp),dimension(:),allocatable :: xabsc,weigh
  func_der_r=0.d0
  func_1_r=0.d0
  count_localx=0
  func_der_theta=0.d0
  func_1_theta=0.d0
  func_2_theta=0.d0
  count_localy=0
  func_der_phi=0.d0
  func_1_phi=0.d0
  count_localz=0

  allocate(xabsc(num_pts))
  allocate(weigh(num_pts))
  call gauleg(num_pts,xabsc,weigh)

  pm=(a+b)/2.d0
  pl=(b-a)/2.d0

  ! R-integral
  do ii=1,num_pts
     x=pm+pl*xabsc(ii)
     rad=radius_tab(x,element_local)
     count_localx=count_localx+1
     func_der_r=func_der_r+weigh(ii)*der_r(x,k_ri)*der_r(x,k_rj)/(an_r*an_r)*an_r*rad**2
     func_1_r=func_1_r+weigh(ii)*bas_r(x,k_ri)*bas_r(x,k_rj)*an_r
  end do
  func_der_r=func_der_r*pl
  func_1_r=func_1_r*pl
  !---------
  !deallocate(xabsc,weigh)
  !num_pts=4
  !allocate(xabsc(num_pts))
  !allocate(weigh(num_pts))
  !call gauleg(num_pts,xabsc,weigh)
  !---------
  ! Theta-integral
  do ii=1,num_pts
     y=pm+pl*xabsc(ii)
     count_localy=count_localy+1
     sin_theta=sin(angle_theta_tab(y,element_local))
     func_der_theta=func_der_theta+weigh(ii)*der_r(y,k_thetai)*der_r(y,k_thetaj)/(an_theta*an_theta)*an_theta*sin_theta
     func_1_theta=func_1_theta+weigh(ii)*bas_r(y,k_thetai)*bas_r(y,k_thetaj)*an_theta*sin_theta
     func_2_theta=func_2_theta+weigh(ii)*bas_r(y,k_thetai)*bas_r(y,k_thetaj)*an_theta/sin_theta
  end do
  func_der_theta=func_der_theta*pl
  func_1_theta=func_1_theta*pl
  func_2_theta=func_2_theta*pl
  !----------
  !deallocate(xabsc,weigh)
  !num_pts=30
  !allocate(xabsc(num_pts))
  !allocate(weigh(num_pts))
  !call gauleg(num_pts,xabsc,weigh)
  !----------
  ! Phi- integral
  do ii=1,num_pts
     count_localz=count_localz+1
     zeta=pm+pl*xabsc(ii)
     func_der_phi=func_der_phi+weigh(ii)*der_r(zeta,k_phii)*der_r(zeta,k_phij)/(an_phi*an_phi)*an_phi
     func_1_phi=func_1_phi+weigh(ii)*bas_r(zeta,k_phii)*bas_r(zeta,k_phij)*an_phi
  end do
  func_der_phi=func_der_phi*pl
  func_1_phi=func_1_phi*pl
  ! Assembling
  separate_integral=func_der_r*func_1_theta*func_1_phi+func_1_r*func_der_theta*func_1_phi+func_der_phi*func_1_r*func_2_theta

  deallocate(xabsc)
  deallocate(weigh)
end function separate_integral


subroutine sing_elim_theta(theta_in,angle,k_thetai,k_thetaj,k_phii,k_phij,sing_elim_choice_phi)
  !*COMM* Subroutine that chooses the function that is one at the singularity
  use nrtype , only : dbl,Pi,i4b
  use Brep,only : epsilon_custom
  !use potential_data
  implicit none
  character(10),intent(out) :: sing_elim_choice_phi
  real(kind=dbl),intent(in) :: theta_in,angle
  integer(kind=i4b), intent(in) :: k_thetai,k_thetaj,k_phii,k_phij
  integer(kind=i4b) :: sum_kr, sum_ktheta
  if ((k_thetai.eq.k_thetaj)) then
     sum_ktheta = k_thetai+k_thetaj
     if ((theta_in.le.epsilon_custom)) then
        if ((sum_ktheta.eq.2).and.((k_phii.eq.2).or.(k_phii.eq.3).or.(k_phij.eq.2).or.(k_phij.eq.3))) then
           sing_elim_choice_phi='func_sing_der'
        else if ((sum_ktheta.eq.2)) then
           sing_elim_choice_phi='func_sing'
        else
           sing_elim_choice_phi='func'
        end if
     else if ((theta_in.le.(Pi-epsilon_custom))) then
        if ((sum_ktheta.eq.4).and.((k_phii.eq.2).or.(k_phii.eq.3).or.(k_phij.eq.2).or.(k_phij.eq.3))) then
           sing_elim_choice_phi='func_sing_der'
        else if ((sum_ktheta.eq.4)) then
           sing_elim_choice_phi='func_sing'
        else
           sing_elim_choice_phi='func'
        end if
     end if
  else
     sing_elim_choice_phi='func'
  end if


end subroutine sing_elim_theta
