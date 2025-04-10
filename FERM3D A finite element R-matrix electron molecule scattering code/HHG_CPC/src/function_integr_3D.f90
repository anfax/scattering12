!!!*************************************************************
! 文件/File: function_integr_3D.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: function_integr_3D.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


module function_integr_tab
  use Calc_func
  use control
  use open_information
  use nrtype,only : Pi,nep_e,order,dbl
  use over, only : overlap_c
  use potential_interp
  use V_setup, only :x_nuc=>x,y_nuc=>y,z_nuc=>z,nAtom,bscoef_V,bscoef_density,xknot,yknot,zknot,&
       & Latter_cutoff,z_cutoff,bscoef_grad,bscoef_nabla
  use bspline
  use gauss3dint
  use bsp_info, only : nE
  implicit none
  INTEGER, PARAMETER :: I4 = SELECTED_INT_KIND(9)
  integer(kind=i4),parameter,private :: dp = kind(1.0d0)
  integer(kind=i4),parameter,private :: kxord=order,kyord=order,kzord=order
  real(kind=dp),private ::  basis,xx,yy,zz,real
  real(kind=dp),private,parameter :: dist=2.078699d0,third=.33333333333333333333d0,phase_dist=.5d0*Pi
  real(kind=dp),private,parameter ::alpha_slater= .75d0
  real(kind=dp),private :: cos_phi,sin_phi,radius,rad,sin_theta,x2,x4,val,rho
  real(kind=dp),dimension(:),allocatable :: r_glo_in_tab,theta_glo_in_tab
  real(kind=dp),private :: hara
  real(kind=dp),private :: density_cuberoot,chiara
  real(kind=dp), private,parameter :: three_Pi_square_cuberoot=3.093667726280135930969d0,two_over_Pi=.63661977236758134307d0
  real(kind=dp),private :: eta,k_Fermi,k_hara,func_eta,plgndr
  external radius,hara,plgndr
  character*15 :: inter
  !---------------------------
  ! The following functions are integrands for the 3D integrals for all
  ! different partitions and matrices, optimized to store a lot of quantities
  ! and not recalculate them every time
  !---------------------------
contains
  !-------------



  real(kind=dp) function func_closed_closed_overlap_tab(x,y,zeta)
    implicit none
    real(kind=dp),intent(in) :: x,y,zeta
    real(kind=dp) :: val_exc,val_V
    rad=radius_tab(x,element)
    num=num+1
    count_pts=count_pts+1
    if (count_func.eq.1) then
       val_exc=interp_density(element,num)

       density_cuberoot=val_exc
       k_Fermi=three_Pi_square_cuberoot*(density_cuberoot)
       k_hara=sqrt(k_Fermi*k_Fermi+2*Energy_last_bound)
       eta=k_hara/k_Fermi
       func_eta=(1.d0/(2.d0*eta)-(1.d0+eta**2)/(4.d0*eta**2)*log((1.0d0+eta)/(eta-1.d0)))/(k_Fermi*k_hara)
       chiara= - two_over_Pi*k_Fermi*func_eta
       !*******
       ! Modified potential
       !chiara=0.d0
       !*******
       val=chiara
       val_array(count_pts)=val
!write(90,30)val,rad,angle_phi(y),angle_theta_tab(zeta,element)
30 format(4e16.8)
    else
       val=val_array(count_pts)
    end if
    func_closed_closed_overlap_tab=(1-val)*mem_tab(count_func,count_pts)*radius_tab(x,element)*sin(angle_theta_tab(zeta,element))*radius_tab(x,element)*an_r*an_theta*an_phi
  end function func_closed_closed_overlap_tab
  !-------------


  real(kind=dp) function func_open_closed_overlap_tab(x,y,zeta)
    implicit none
    real(kind=dp),intent(in) :: x,y,zeta
    real(kind=dp) :: val_exc,val_V
    num=num+1
    count_pts=count_pts+1
    val_exc=interp_density(element,num)

    density_cuberoot=val_exc
    k_Fermi=three_Pi_square_cuberoot*(density_cuberoot)
    k_hara=sqrt(k_Fermi*k_Fermi+2*Energy_last_bound)
    eta=k_hara/k_Fermi
    func_eta=(1.d0/(2.d0*eta)-(1.d0+eta**2)/(4.d0*eta**2)*log((1.0d0+eta)/(eta-1.d0)))/(k_Fermi*k_hara)
    chiara= - two_over_Pi*k_Fermi*func_eta

    !******
    ! Modified potential


    !chiara=0.d0
    !******

    val=chiara

    func_open_closed_overlap_tab=(1-val)*((bas_tab(k_ri,countx)*bas_tab(k_rj,countx)*bas_tab(k_phii,county)*bas_tab(k_phij,county)*&
         &bas_tab(k_thetai,countz)*bas_tab(k_thetaj,countz))*radius_tab(x,element)*sin(angle_theta_tab(zeta,element))*radius_tab(x,element))*an_r*an_theta*an_phi
  end function func_open_closed_overlap_tab
  !-------------



  !--------------
  real(kind=dp) function func_closed_closed_potential_tab(x,y,zeta)
    implicit none
    real(kind=dp),intent(in) :: x,y,zeta
    real(kind=dp) :: val_exc,val_V
    rad=radius_tab(x,element)
    num=num+1
    count_pts=count_pts+1
    if (count_func.eq.1) then
       val_V=interp_V(element,num)
       val_exc=interp_density(element,num)

       density_cuberoot=val_exc
       k_Fermi=three_Pi_square_cuberoot*(density_cuberoot)
       k_hara=sqrt(k_Fermi*k_Fermi+2*Energy_last_bound)
       eta=k_hara/k_Fermi
       func_eta=.5d0+((1.0d0 -eta*eta)/(4.d0*eta))*log(abs((1.d0+eta)/(1.d0-eta)))
       chiara= - two_over_Pi*k_Fermi*func_eta

       val=val_V+chiara
       !*******
       ! Modified potential
       !val=-1.d0/(sqrt(rad**2+dist**2-2*rad*dist*cos(angle_theta_tab(zeta,element))))
       !*******
       !write(1001,*)val,rad,angle_theta_tab(zeta,element)
       !write(1002,*)val,rad,angle_phi(y)
       !---------------------------------------------------------------------------
       ! Get rid of it
           !xx=radius(x)*sin(angle_theta_tab(zeta,element))*cos(angle_phi(y))
           !yy=radius(x)*sin(angle_theta_tab(zeta,element))*sin(angle_phi(y))
           !zz=radius(x)*cos(angle_theta_tab(zeta,element))
       !if ((abs(xx).le. 0.025d0).and.(abs(yy).le. 0.025d0)) then
       !write(1001,900)val_V,zz,xx,yy
       !write(1002,900)chiara,zz,xx,yy
       !write(1003,900)val-chiara-val_V,zz,xx,yy
       !end if
       !900 format(4e16.8,2I4)
       !---------------------------------------------------------------------------
       !---------
       if ((.not.DFT)) then
          val=val-(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)*(1- exp(-(rad/r_cutoff)**6))
       else if (abs(density_cuberoot).lt.1.d-12) then
          val=val-(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)
          !---------
       else if (rad.lt.r_cutoff) then
          inter='outside'
          val=val+DFT_correl(radius_tab(x,element),density_cuberoot,r_cutoff,zeta,num,element,y,inter)
       else
          val=val-(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)
       end if

       val_array(count_pts)=val
!if (abs(xx).lt.0.3d0)  write(91,30)val,xx,yy,zz
!if (abs(yy).lt.0.3d0)  write(92,30)val,xx,yy,zz

!write(91,30)val,rad,angle_phi(y),angle_theta_tab(zeta,element)
30 format(4e16.8)

    else
       val=val_array(count_pts)
    end if
    func_closed_closed_potential_tab=(val)*mem_tab(count_func,count_pts)*rad*sin(angle_theta_tab(zeta,element))*rad*an_r*an_theta*an_phi

  end function func_closed_closed_potential_tab
  !--------------

  real(kind=dp) function func_open_closed_potential_tab(x,y,zeta)
    implicit none
    real(kind=dp),intent(in) :: x,y,zeta
    real(kind=dp) :: val_exc,val_V
    rad=radius_tab(x,element)
    num=num+1
    val_V=interp_V(element,num)
    val_exc=interp_density(element,num)
    density_cuberoot=val_exc
    k_Fermi=three_Pi_square_cuberoot*(density_cuberoot)
    k_hara=sqrt(k_Fermi*k_Fermi+2*Energy_last_bound)
    eta=k_hara/k_Fermi
    func_eta=.5d0+((1.0d0 -eta*eta)/(4.d0*eta))*log(abs((1.d0+eta)/(1.d0-eta)))
    chiara= - two_over_Pi*k_Fermi*func_eta
    !val=val_V+chiara -(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)*(1- exp(-(rad/r_cutoff)**6))
    !-----------
    if ((.not.DFT)) then
       val=val_V+chiara -(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)*(1- exp(-(rad/r_cutoff)**6))
    else if (abs(density_cuberoot).lt.1.d-12) then
       val=val_V+chiara -(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)
    else 
       val=val_V+chiara -(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)
    end if
    !-----------
    !else
    !inter='outside'
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! DEBUG 
    !val=val_V+chiara+DFT_correl(radius(x),density_cuberoot,r_cutoff,zeta,num,element,y,inter)
    !val=val_V+chiara -(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)*(1- exp(-(rad/r_cutoff)**6))
    !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !end if
    !********
    ! Modified potential
    !val=-1.d0/(sqrt(rad**2+dist**2-2*rad*dist*cos(angle_theta_tab(zeta,element))))
    !******** 
    func_open_closed_potential_tab=val*((bas_tab(k_ri,countx)*bas_tab(k_rj,countx)*bas_tab(k_phii,county)*bas_tab(k_phij,county)*bas_tab(k_thetai,countz)*bas_tab(k_thetaj,countz))*rad*sin(angle_theta_tab(zeta,element))*rad)*an_r*an_theta*an_phi
  end function func_open_closed_potential_tab
  !-------------

  real(kind=dp) function func_sing_1_tab(x,y,zeta)
    implicit none
    real(kind=dp),intent(in) :: x,y,zeta
    real(kind=dp) :: val_exc,val_V
    num=num+1
    rad=radius_tab(x,element)
    sin_theta=sin(angle_theta_tab(zeta,element))
    count_pts=count_pts+1
    if (count_func.eq.1) then
       val_V=interp_V_sing(ind_sing,num)
       val_exc=interp_density_sing(ind_sing,num)
       density_cuberoot=val_exc
       k_Fermi=three_Pi_square_cuberoot*(density_cuberoot)
       k_hara=sqrt(k_Fermi*k_Fermi+2*Energy_last_bound)
       eta=k_hara/k_Fermi
       func_eta=.5d0+((1.0d0 -eta*eta)/(4.d0*eta))*log(abs((1.d0+eta)/(1.d0 -eta)))
       chiara= - two_over_Pi*k_Fermi*func_eta
       !val=val_V+chiara-(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)*(1-exp(-(rad/r_cutoff)**6))
       if ((.not.DFT)) then
          val=val_V+chiara-(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)*(1-exp(-(rad/r_cutoff)**6))
       else if (abs(density_cuberoot).lt.1.d-12) then
          val=val_V+chiara-(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)
       else if (rad.lt.r_cutoff) then
          !---------
          inter='sing'
          ! DEBUG
          val=val_V+chiara+DFT_correl(radius_tab(x,element),density_cuberoot,r_cutoff,zeta,num,element,y,inter)
       else
          val=val_V+chiara-(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,element))))/(rad**4)
       end if
       !******
       !val=-1.d0/(sqrt(rad**2+dist**2-2*rad*dist*cos(angle_theta_tab(zeta,element))))
       !******
       !write(3001,*)val,rad,angle_theta_tab(zeta,element)
       !write(3002,*)val,rad,angle_phi(y)
       !---------------------------------------------------------------------------
       ! Get rid of it
       !    xx=radius(x)*sin(angle_theta_tab(zeta,element))*cos(angle_phi(y))
       !    yy=radius(x)*sin(angle_theta_tab(zeta,element))*sin(angle_phi(y))
       !    zz=radius(x)*cos(angle_theta_tab(zeta,element))
       !if ((abs(xx).le. 0.025d0).and.(abs(yy).le. 0.025d0)) then
       !write(1001,900)val_V,zz,xx,yy
       !write(1002,900)chiara,zz,xx,yy
       !write(1003,900)val-chiara-val_V,zz,xx,yy
       !end if
       !900 format(4e16.8,2I4)
       !---------------------------------------------------------------------------
       val_array_sing(count_pts)=val
    else
       val=val_array_sing(count_pts)
    end if
    func_sing_1_tab=val*((bas_tab_sing(k_ri,countx)*bas_tab_sing(k_rj,countx)*bas_tab_sing(k_phii,county)*bas_tab_sing(k_phij,county)*bas_tab_sing(k_thetai,countz)*bas_tab_sing(k_thetaj,countz)))*rad*sin_theta*rad*an_r*an_theta*an_phi
  end function func_sing_1_tab
  !-------------


  real(kind=dp) function func_sing_0_tab(x,y,zeta)
    implicit none
    real(kind=dp),intent(in) :: x,y,zeta
    real(kind=dp) :: val_exc,val_V
    num=num+1
    rad=radius_tab(x,element)
    sin_theta=sin(angle_theta_tab(zeta,element))
    count_pts=count_pts+1
    if (count_func.eq.1) then
       val_V=interp_V_sing(ind_sing,num)
       val_exc=interp_density_sing(ind_sing,num)
       density_cuberoot=val_exc
       k_Fermi=three_Pi_square_cuberoot*(density_cuberoot)
       k_hara=sqrt(k_Fermi*k_Fermi+2*Energy_last_bound)
       eta=k_hara/k_Fermi
       func_eta=(1.d0/(2.d0*eta)-(1.d0+eta**2)/(4.d0*eta**2)*log((1.0d0+eta)/(eta-1.d0)))/(k_Fermi*k_hara)
       chiara= - two_over_Pi*k_Fermi*func_eta
       !*****
       ! Modified potential

       !chiara=0.d0
       !*****
       !write(4001,*)val,rad,angle_theta_tab(zeta,element)
       !write(4002,*)val,rad,angle_phi(y)
       val=chiara
       val_array_sing(count_pts)=val
    else
       val=val_array_sing(count_pts)
    end if

    func_sing_0_tab=(1-val)*((bas_tab_sing(k_ri,countx)*bas_tab_sing(k_rj,countx)*bas_tab_sing(k_phii,county)*bas_tab_sing(k_phij,county)*bas_tab_sing(k_thetai,countz)*bas_tab_sing(k_thetaj,countz)))*rad*sin_theta*rad*an_r*an_theta*an_phi
  end function func_sing_0_tab
  !-------------


  real(kind=dp) function func_potential_empty_tab(x,y,zeta)
    implicit none
    real(kind=dp),intent(in) :: x,y,zeta
    real(kind=dp) :: val_exc,weater_potential,val_grad,val_nabla
    external :: weater_potential
    num=num+1
    xx=radius(x)*sin(angle_theta(zeta))*cos(angle_phi(y))
    yy=radius(x)*sin(angle_theta(zeta))*sin(angle_phi(y))
    zz=radius(x)*cos(angle_theta(zeta))
    rho=sqrt(xx**2+yy**2)
    val=dbs3vl(xx,yy,zz,kxord,kyord,kzord,xknot,yknot,zknot,nnxdata,nnydata,nnzdata,bscoef_V)
    do i=1,nAtom
       val =val/sqrt((xx-x_nuc(i))**2+(yy-y_nuc(i))**2+(zz-z_nuc(i))**2)
    end do
    val_exc=dbs3vl(xx,yy,zz,kxord,kyord,kzord,xknot,yknot,zknot,nnxdata,nnydata,nnzdata,bscoef_density)
    interp_V(element,num)=val
    interp_density(element,num)=val_exc
    func_potential_empty_tab=0.0d0
    if (DFT) then
       val_grad=dbs3vl(xx,yy,zz,kxord,kyord,kzord,xknot,yknot,zknot,nnxdata,nnydata,nnzdata,bscoef_grad)
       val_nabla=dbs3vl(xx,yy,zz,kxord,kyord,kzord,xknot,yknot,zknot,nnxdata,nnydata,nnzdata,bscoef_nabla)
       interp_grad(element,num)=val_grad
       interp_nabla(element,num)=val_nabla
    end if

  end function func_potential_empty_tab
  !-------------


  real(kind=dp) function func_potential_empty_sing_tab(x,y,zeta)
    implicit none
    real(kind=dp),intent(in) :: x,y,zeta
    real(kind=dp) :: val_exc,weater_potential,val_grad,val_nabla
    external :: weater_potential
    num=num+1
    xx=radius(x)*sin(angle_theta(zeta))*cos(angle_phi(y))
    yy=radius(x)*sin(angle_theta(zeta))*sin(angle_phi(y))
    zz=radius(x)*cos(angle_theta(zeta))
    rho=sqrt(xx**2+yy**2) 
    val=dbs3vl(xx,yy,zz,kxord,kyord,kzord,xknot,yknot,zknot,nnxdata,nnydata,nnzdata,bscoef_V)
    do i=1,nAtom
       val =val/sqrt((xx-x_nuc(i))**2+(yy-y_nuc(i))**2+(zz-z_nuc(i))**2)
    end do
    val_exc=dbs3vl(xx,yy,zz,kxord,kyord,kzord,xknot,yknot,zknot,nnxdata,nnydata,nnzdata,bscoef_density)
    interp_V_sing(ind_sing,num)=val
    interp_density_sing(ind_sing,num)=val_exc
    func_potential_empty_sing_tab=0.0d0
    if (DFT) then
       val_grad=dbs3vl(xx,yy,zz,kxord,kyord,kzord,xknot,yknot,zknot,nnxdata,nnydata,nnzdata,bscoef_grad)
       val_nabla=dbs3vl(xx,yy,zz,kxord,kyord,kzord,xknot,yknot,zknot,nnxdata,nnydata,nnzdata,bscoef_nabla)
       interp_grad_sing(ind_sing,num)=val_grad
       interp_nabla_sing(ind_sing,num)=val_nabla
    end if
  end function func_potential_empty_sing_tab
  !-------------

  real(kind=dp) function func_mem_tab(x,y,zeta)
    implicit none
    real(kind=dp),intent(in) :: x,y,zeta
    count_pts=count_pts+1
    mem_tab(count_func,count_pts)=(bas_tab(k_ri,countx)*bas_tab(k_rj,countx)*bas_tab(k_phii,county)*bas_tab(k_phij,county)*bas_tab(k_thetai,countz)*bas_tab(k_thetaj,countz)) 
    func_mem_tab=0.d0
  end function func_mem_tab

  !--------------
  real(kind=dp) function func_tab(x,y,zeta)
    implicit none
    real(kind=dp),intent(in) :: x,y,zeta
    integer(kind=i4) k_x
    num=num+1
    do k_x=1,4
       bas_tab(k_x,countx)=bas_r(x,k_x)
       der_tab(k_x,countx)=der_r(x,k_x)
    end do
    func_tab=0.0d0
  end function func_tab

  real(kind=dp) function func_tab_sing(x,y,zeta)
    implicit none
    real(kind=dp),intent(in) :: x,y,zeta
    integer(kind=i4) k_x
    num=num+1
    do k_x=1,4
       bas_tab_sing(k_x,countx)=bas_r(x,k_x)
       der_tab_sing(k_x,countx)=der_r(x,k_x)
    end do
    func_tab_sing=0.0d0
  end function func_tab_sing


  real(kind=dp) function radius_tab(x,element)
    use nrtype,only : dp
    use Brep
    use Calc_func,only : an_r
    use control
    use open_information
    implicit none
    integer(kind=i4) element
    real(kind=dp),intent(in) :: x
    real(kind=dp) :: r_glo_in
    external r_glo_in
    radius_tab=r_glo_in_tab(element)+x*an_r
  end function radius_tab


  real(kind=dp) function angle_theta_tab(y,element)
    !use nrtype,only : dbl
    use Brep
    use Calc_func,only :an_theta
    use control
    use open_information
    implicit none
    integer(kind=i4) element
    real(kind=dp),intent(in) :: y
    real(kind=dp) :: theta_glo_in
    external theta_glo_in
    angle_theta_tab=theta_glo_in_tab(element)+y*an_theta
  end function angle_theta_tab

  subroutine r_glo_in_t
    use Brep
    !        use nrtype, only : dbl
    implicit none
    integer(kind=i4) :: l
    integer(kind=i4):: element
    integer(kind=i4),parameter :: first=1
    real(kind=dp) :: r_glo_in_temp
    allocate(r_glo_in_tab(n))
    ! This subroutine calculates the initial radius value for each element/sector
    do element=1,n  
       r_glo_in_temp= r_glo(nodes(element,first))
       do l=1,8
          if (r_glo(nodes(element,l)).lt.r_glo_in_temp) then
             r_glo_in_temp=r_glo(nodes(element,l))
          end if
       end do
       r_glo_in_tab(element)=r_glo_in_temp
    end do
  end subroutine r_glo_in_t


  subroutine theta_glo_in_t
    use Brep
    !        use nrtype, only : dbl
    implicit none
    integer(kind=i4) :: l
    integer(kind=i4):: element
    integer(kind=i4),parameter :: first=1
    real(kind=dp) :: theta_glo_in_temp
    allocate(theta_glo_in_tab(n))
    ! This subroutine calculates the initial radius value for each element/sector
    do element=1,n
       theta_glo_in_temp= theta_glo(nodes(element,first))
       do l=1,8
          if (theta_glo(nodes(element,l)).lt.theta_glo_in_temp) then
             theta_glo_in_temp=theta_glo(nodes(element,l))
          end if
       end do
       theta_glo_in_tab(element)=theta_glo_in_temp
    end do
  end subroutine theta_glo_in_t


  function distance_phi(r,rr,t,tt,f,ff)
    use nrtype
    implicit none
    real (kind=dbl),intent (in) :: r,rr,t,tt,f,ff 
    real (kind=dbl) distance_phi
    distance_phi=Sqrt(r**2 + rr**2 - 2.d0*r*rr*(Cos(t)*Cos(tt) + Cos(f - ff)*Sin(t)*Sin(tt)))

  end function distance_phi

  real(kind=dp) function DFT_correl(rad,dens,rcutoff,zeta,num,elem,y,inter)
    real(kind=dp),parameter :: three_quarter_cuberoot=1.330670039491469d0,a=0.04918d0,b=0.132d0,c=0.2533d0,d=0.349d0,C_f=2.871234000188191d0
    real(kind=dp),intent(in) :: rad,dens,rcutoff,zeta,y
    integer, intent(in) :: num,elem
    character*15,intent(in) :: inter
    integer :: na
    real(kind=dp) :: rs,F,G,Fp,Gp,Fd,Gd,grad,nabla,V_pol
    !if (rad.lt.rcutoff) then
       xx=rad*sin(angle_theta_tab(zeta,elem))*cos(angle_phi(y))
       yy=rad*sin(angle_theta_tab(zeta,elem))*sin(angle_phi(y))
       zz=rad*cos(angle_theta_tab(zeta,elem))
       F = 1.d0/(1.d0+d/dens)
       G = F*dens**(-5)*exp(-c/dens)
       Fp = d/(3.d0*(1.d0 + d/dens)**2*dens**4)
       Gp =    (c*d + (c - 4.d0*d)*dens- 5.d0*dens**2)/(3.d0*exp(c/dens)*(d +dens)**2*dens**8)
       Fd = (-2*d*(d + 2.d0*dens))/(9.d0*(d + dens)**3*dens**5)
       Gd =     (c**2*d**2 + 2.d0*c*(c - 6.d0*d)*d*dens+ (c**2 - 26.d0*c*d + 28.d0*d**2)*dens**2 &
            &-  14.d0*c*dens**3 + 66.d0*d*dens**3 + 40.d0*dens**4)/(9.d0*exp(c/dens)*(d +dens)**3*dens**12)
       !-------------------------------

       if (inter.eq.'sing') then
          grad=(interp_grad_sing(ind_sing,num))
          nabla=(interp_nabla_sing(ind_sing,num))
       else
          grad=(interp_grad(elem,num))
          nabla=(interp_nabla(elem,num))
       end if
       do na=1,nAtom
          nabla=nabla / ((xx-x_nuc(na))**2+(yy-y_nuc(na))**2+(zz-z_nuc(na))**2)
       end do
       do na=1,nAtom
          grad=grad / ((xx-x_nuc(na))**2+(yy-y_nuc(na))**2+(zz-z_nuc(na))**2)
       end do

       !grad=grad**2

       !-------------------------------
       ! Density Functional from BLYP
       DFT_correl=-a * (Fp*dens**3+F)- a*b* C_f*dens**5 *(Gp*dens**3+8.d0/3.d0*G)-a*b/4*(Gd*dens**3*&
            & grad**2+Gp*(3.d0*grad**2+2.d0*dens**3*nabla)+4.d0*G*nabla)-   a*b/72*(3*Gd* dens**3*grad**2+ &
            & Gp*(5.d0*grad**2+6*dens**3*nabla)+4*G*nabla)
       DFT_correl=2*DFT_correl
       if (DFT_correl.gt.0.d0) then
          DFT_correl=0.d0
       end if
    V_pol=-(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,elem))))/(rad**4)
    if (DFT_correl.lt.V_pol) DFT_correl=V_pol

       !-------------------------------
    !else
    !   DFT_correl=-(polar+polar_2*plgndr(2,0,cos(angle_theta_tab(zeta,elem))))/(rad**4)
    !end if
    !XXXXXXXXXXXXXX
    !DEBUG
    if (angle_theta_tab(zeta,elem).lt..1d0) then
       !write(4001,900)DFT_correl,xx,yy,zz
       !write(196,900)G,Gp,Gd,rad
       write(197,*)dens,rad
       !write(194,900)F,Fp,Fd,rad
       write(200,*)DFT_correl,rad,angle_phi(y)
       !write(199,*)V_pol,rad,angle_phi(y)
       write(201,*)grad,rad
       write(202,*)nabla,rad,inter
       !DFT_correl=0.d0
    end if
    !XXXXXXXXXXXXXX
900 format(4e14.6)

  end function DFT_correl


end module function_integr_tab


