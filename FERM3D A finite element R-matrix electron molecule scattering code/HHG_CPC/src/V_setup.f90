!!!*************************************************************
! 文件/File: V_setup.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: V_setup.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


module V_setup
  use nrtype,only : dbl,i4b
  real(kind=dbl),dimension(:),allocatable :: Charge,X,Y,Z
  real(kind=dbl),dimension(:),allocatable :: position,angle,angle_position_phi
  real(kind=dbl),dimension(:,:,:),allocatable :: bscoef_V,bscoef_density,bscoef_grad,bscoef_nabla
  real(kind=dbl),dimension(:),allocatable :: xknot,yknot,zknot
  real(kind=dbl),dimension(:),allocatable :: xx,yy,zz
  real(kind=dbl) :: Latter_cutoff=5.0d0,z_cutoff=5.5d0,molec_lin
  integer(kind=i4b) :: nAtom,n1,n2,n3
contains

  subroutine V_conversion(n,nx)
    use Brep,only : r_glo_avg,phi_glo_avg,theta_glo_avg,epsilon_custom,R0,shift
    use nrtype, only : Pi,dbl,order,nep_e,i4b
    !    use bsp_info 
    use bspline , only : dbs3in,dbs3vl
    use control , only : DFT,restore
    implicit none
    interface 
       subroutine dbs3ines(xdata,ydata,zdata,fdata,nxdata,nydata,nzdata,nxknot,nyknot,nzknot,bscoef,xknot,yknot,zknot)
         integer nxdata,nydata,nzdata,nxknot,nyknot,nzknot
         real(kind=8), dimension(nxdata),intent(in) :: xdata
         real(kind=8), dimension(nydata),intent(in) :: ydata
         real(kind=8), dimension(nzdata),intent(in) :: zdata
         real(kind=8), dimension(nxdata,nydata,nzdata), intent(in) :: fdata
         real(kind=8), dimension(nxdata,nydata,nzdata), intent(out) :: bscoef
         real(kind=8), dimension(nxknot),intent(out) :: xknot
         real(kind=8), dimension(nyknot),intent(out) :: yknot
         real(kind=8), dimension(nzknot),intent(out) :: zknot
       end subroutine dbs3ines
    end interface
    integer(kind=i4b),intent(in) :: n,nx
    integer(kind=i4b),parameter  :: ny = 15
    real(kind=dbl),dimension(:,:,:),allocatable :: V,density
!    real(kind=dbl),dimension(:),allocatable :: xx,yy,zz
    real(kind=dbl),dimension(:,:),allocatable :: V_t
    integer(kind=i4b), dimension(:),allocatable ::Atomic_num
    integer(kind=i4b) i1,i2,i3,element,i,na
    integer(kind=i4b) aliment
    integer(kind=i4b) nxknot,nyknot,nzknot
    real(kind=dbl) :: distance_temp,distance_x,x0,y0,z0,x1,x2,x3,y1,y2,y3,z1,z2,z3,val,x_coord,y_coord,z_coord,alpha
    real(kind=dbl),dimension(:),allocatable :: x_coordinate,y_coordinate,z_coordinate
    character*40 :: title,title2
    real(kind=dbl),parameter :: third = .3333333333333333333333333d0 
    integer(kind=i4b),parameter ::m=0 
    real(kind=dbl) :: raggio,plgndr,cos_theta,rho,delta,theta,temp,shifty,shiftx,beta
    real(kind=dbl),parameter :: pre_V=.5d0,pre_density=1.0d-7,lower_radius_limit=.5d0,epsilon_dens=3.0d-7
    real(kind=dbl),dimension(:), allocatable :: theta_d,raggio_d,rho_d,z_d 
    integer(kind=i4b) :: test,l,k,ios,ios_rot,ios_DFT,ios_x,ios_restore
    integer(kind=i4b), parameter :: num_theta_sectors=100,nw=12
    logical :: rotation=.FALSE.

    !****************************************************
    ! Subroutine that interpolates electrostatis potential and electron density
    ! 
    ! Variables
    ! Order: order of the splines in 'interpolation_imsl'
    ! nx-y-zknot : dimension of the knot vector for the spline interp.
    ! 
    !****************************************************

    write(6,*)'order interpolant splines=',order
    !Reading the potential from the cube file

    !$write(32,*)VV
    open(unit=nx,file='pot.dat',status='old',action='read')
    open(unit=ny,file='dens.dat',status='old',action='read')
    read(nx,9000)title
    read(nx,9000)title2
9000 format(1A40)
    read(nx,*)nAtom,x0,y0,z0
    write(6,*)'number atoms, first grid point',nAtom,x0,y0,z0

    !*COMM***
    ! Shift for elimination of dipole moment for eteronuclear diatomics
    ! The atom that presumably brings the most of the charge (for ions) is being
    ! put on the origin
    !************
    z0=z0-shift

    open(unit=nw,file='input_control.dat',status='old',action='read')
    do i=1,12
    read(nw,*)
    end do
    read(nw,*,IOSTAT=ios)shifty
    write(6,*)'IOSTAT=',ios,shifty
    read(nw,*,IOSTAT=ios_rot)rotation
    write(6,*)'IOSTAT=',ios_rot,rotation
    read(nw,*,IOSTAT=ios_DFT)DFT
    write(6,*)'IOSTAT=',ios_DFT,DFT

!-------------------------------------
    read(nw,*,IOSTAT=ios_x)shiftx
    write(6,*)'IOSTAT=',ios_x,shiftx
    read(nw,*,IOSTAT=ios_restore)restore
    write(6,*)'IOSTAT=',ios_restore,restore
    close(nw)   ! Close input_control.dat file
    !************
    !*COMM***
    ! Shift of coordinates on the x axis
    if (ios_x.eq.0) then
       x0=x0-shiftx
       write(6,*)'molecule shifted on x of',shiftx
    end if
    !************


!-------------------------------------

    !*COMM***
    ! Shift of coordinates on the y axis
    if (ios.eq.0) then
       y0=y0-shifty
       write(6,*)'molecule shifted on y of',shifty
    end if
    !************
    allocate(Atomic_num(nAtom))
    allocate(Charge(nAtom))
    allocate(position(nAtom))
    allocate(angle(nAtom))
    allocate(angle_position_phi(nAtom))
    allocate(x(nAtom))
    allocate(y(nAtom))
    allocate(z(nAtom))
    read(nx,*)n1,x1,y1,z1
    read(nx,*)n2,x2,y2,z2
    read(nx,*)n3,x3,y3,z3
       write(6,*)'atomic number,charge,coordinates '
    do i=1,nAtom
       read(nx,*)Atomic_Num(i),Charge(i),x(i),y(i),z(i)
       !*************
       z(i)=z(i)-shift
       !---------------------
       ! Shift of coordinates on the y axis
       if (ios.eq.0) y(i)=y(i)-shifty
       !---------------------
       !---------------------
       ! Shift of coordinates on the x axis
       if (ios_x.eq.0) x(i)=x(i)-shiftx
       !---------------------
       ! Pseudo rotation of the axes
       if ((ios_rot.eq.0).and.rotation) then
          temp=y(i)
          y(i)=z(i)
          z(i)=temp
       end if
       !---------------------
       !*************
       write(6,901)Atomic_Num(i),Charge(i),x(i),y(i),z(i)
    end do
    !-------------
    !    read(ny,*)title
    !    read(ny,*)title2
    !    read(ny,*)nAtom,x0,y0,z0
    !    read(ny,*)n1,x1,y1,z1
    !    read(ny,*)n2,x2,y2,z2
    !    read(ny,*)n3,x3,y3,z3
    !    do i=1,nAtom
    !       read(ny,*)Atomic_Num(i),Charge(i),x(i),y(i),z(i)
    !       write(6,*)Atomic_Num(i),Charge(i),x(i),y(i),z(i)
    !    end do
    do i=1,6+nAtom
       read(ny,*)
    end do
    !-------------


    !-------------
    !*COMM* get positions of nuclei in spherical coordinates to subtract
    !*COMM* singularities
       write(6,*)"position of atoms in the grid: r, theta, phi"
    do i=1,nAtom
       position(i)=Sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
       if (position(i).gt. 0.0d0+epsilon_custom) then
          angle(i)=acos(z(i)/position(i))
          if ((angle(i).gt.0.d0+epsilon_custom).and.(angle(i).le.Pi-epsilon_custom)) then
             beta=x(i)/(position(i)*sin(angle(i)))
             if (beta.le.-1.d0) then
                angle_position_phi(i)=Pi
             else if (beta.gt.1.d0) then
                angle_position_phi(i)=0.d0
             else
                angle_position_phi(i) = dacos(beta)
             end if
             if (y(i).lt.0.d0) then
                angle_position_phi(i)=2.d0*Pi-angle_position_phi(i)
             end if
          else
             angle_position_phi(i) = 4.d0 * Pi
          end if
       else
          !*COMM* For atoms in the center of the grid
          angle_position_phi(i) = 4.d0 * Pi
          angle(i) = 2.0d0*Pi
       end if
       write(6,901)i,position(i),angle(i),angle_position_phi(i)
    end do
    !-------------
    allocate(V(n1,n2,n3))
    allocate(density(n1,n2,n3))
    do i1=1,n1
       do i2=1,n2
          read(nx,900)( V(i1,i2,i3),i3=1,n3)
          read(ny,900)( density(i1,i2,i3),i3=1,n3)
       end do
    end do
900 format(6E13.5)
901 format(1I5,6E13.5)


    !*COMM*
    !YYYYYYY WARNING :  Flip sign to the potential 
    ! It comes out of gaussian with the opposite sign
    V=-V
    !YYYYYYY END OF WARNING
    !*COMM*

    allocate(xx(n1))
    allocate(yy(n2))
    allocate(zz(n3))
    ! Interpolation (cubic splines)


    do i1=1,n1
       xx(i1) = x0+(i1-1)*x1
    end do
    do i2=1,n2
       yy(i2) = y0+(i2-1)*y2
    end do
    do i3=1,n3
       zz(i3) = z0+(i3-1)*z3
    end do
    !---------------
    ! Implementation of the Latter tail 1- Spherical surface
    if ((ios_rot.eq.0).and.rotation) then
       allocate(V_t(n2,n3))
       do i1=1,n1
          V_t(:,:)=V(i1,:,:)
          V_t=transpose(V_t)
          V(i1,:,:)=V_t(:,:)
          V_t(:,:)=density(i1,:,:)
          V_t=transpose(V_t)
          density(i1,:,:)=V_t(:,:)
       end do
       deallocate(V_t)
    end if
    do i1=1,n1
       do i2=1,n2
          do i3=1,n3
             raggio=sqrt(xx(i1)**2+yy(i2)**2+zz(i3)**2)
             cos_theta=(zz(i3)/raggio)
             if (abs(yy(i2)).le..1d0) then
                if (abs(xx(i1)).le..1d0) then
                   !if ((abs(zz(i3)).le. 4.d0).and.(abs(zz(i3)).ge..0d0)) then
                end if
             end if
             if (V(i1,i2,i3).gt. 1.d9) then
                write(6,*)'Error in interpolation of electrostatic potential'
                stop
             end if
          end do
       end do
    end do
    !----------------
    ! Implementation of the Latter tail 2- Cylindrical surface
    !       do i1=1,n1
    !          do i2=1,n2
    !             do i3=1,n3
    !                rho=sqrt(xx(i1)**2+yy(i2)**2)
    !                if ((rho.ge.Latter_cutoff).or.(zz(i3).ge.z_cutoff)) then
    !                   raggio=sqrt(xx(i1)**2+yy(i2)**2+zz(i3)**2)
    !                   V(i1,i2,i3)=.5d0*((nep_e**(-(raggio/Latter_cutoff)**6)))
    !                end if
    !             end do
    !          end do
    !       end do
    !----------------
    ! *COMM* Theta-dependent cutoff
    allocate(theta_d(num_theta_sectors))
    allocate(raggio_d(num_theta_sectors))
    allocate(rho_d(num_theta_sectors))
    allocate(z_d(num_theta_sectors))
    delta=Pi/num_theta_sectors
    do k=1,num_theta_sectors
       theta_d(k)=k*delta
       raggio_d(k)=R0
    end do
    do i1=1,n1
       do i2=1,n2
          do i3=1,n3
             if (density(i1,i2,i3).lt.epsilon_dens ) then
                raggio=sqrt(xx(i1)**2+yy(i2)**2+zz(i3)**2)
                theta=acos(zz(i3)/raggio)

                if (raggio.ge.lower_radius_limit) then
                   do k=1,num_theta_sectors
                      if ((theta.le.(theta_d(k))).and.(theta.gt.(theta_d(k)-delta))) then
                         if (raggio.lt.raggio_d(k)) then
                            raggio_d(k)=raggio
                            rho_d(k)=sqrt(xx(i1)**2+yy(i2)**2)
                            z_d(k)=zz(i3) 
                         end if
                      end if
                   end do
                end if
             end if
          end do
       end do
    end do
    ! Plotting of potential and density
    do i1=1,n1
       do i2=1,n2
          do i3=1,n3
             raggio=sqrt(xx(i1)**2+yy(i2)**2+zz(i3)**2)
             theta=acos(zz(i3)/raggio)
             if (raggio.ge.lower_radius_limit) then
                do k=1,num_theta_sectors
                   if ((theta.le.(theta_d(k))).and.(theta.gt.(theta_d(k)-delta))) then
                      if (raggio_d(k).le.raggio) then
                         !V(i1,i2,i3)=pre_V*((nep_e**(-(raggio/raggio_d(k))**6)))
                         !density(i1,i2,i3)=pre_density*((nep_e**(-(raggio/raggio_d(k))**6)))
                         !V(i1,i2,i3)=0.0d0
                         !density(i1,i2,i3)=0.0d0
                         !write(107,*)V(i1,i2,i3),raggio
                         !write(109,*)V(i1,i2,i3),xx(i1),zz(i3)
                      else
                         !write(107,*)V(i1,i2,i3),raggio
                         !write(109,*)V(i1,i2,i3),xx(i1),zz(i3)
                      end if
                   end if
                end do
             end if
          end do
       end do
    end do
    !write(6,*)'pre_dens'

    !----------------
    !*COMM* Interpolation of cubic root of density to interpolate a smoother function

    do i1=1,n1
       do i2=1,n2
          do i3=1,n3
             !                write(37,*)density(i1,i2,i3)

             !********************************
             ! Prevention of error for extraction of cubic root of density
             ! density<0 sometimes in the input (abs(density)<1E-10 in these cases)

             if (density(i1,i2,i3).lt.0.d0) then
                density(i1,i2,i3)=0.d0
             end if
             !********************************
             density(i1,i2,i3)=density(i1,i2,i3)**third
          end do
       end do
    end do

    !----------------

    !----------------
    !*COMM* Interpolation of V*r_i*r_j to avoid cusps in the interpolation
    !write(6,*)'pre_pot'
    do i1=1,n1
       do i2=1,n2
          do i3=1,n3
             do na=1,nAtom
                V(i1,i2,i3) = V(i1,i2,i3) * sqrt((xx(i1)-x(na))**2+(yy(i2)-y(na))**2+(zz(i3)-z(na))**2)
             end do
          end do
       end do
    end do

    !----------------
    !*COMM* Plot of modified potential 
    !do i1=1,n1
    !  do i2=1,n2
    !    do i3=1,n3
    !       if (abs(zz(i3)).le..005) then
    !         if (abs(yy(i2)).le..1) then
    !           if (xx(i1).ge.0.0d0) then
    !             write(3005,*)xx(i1),zz(i2),V(i1,i2,i3)
    !         write(3002,*)V(i1,i2,i3),xx(i1),zz(i3)
    !           end if
    !         end if
    !       end if
    !    end do
    !  end do
    !end do

    !----------------
    nxknot=order+n1
    nyknot=order+n2
    nzknot=order+n3
    allocate(bscoef_V(n1,n2,n3))
    allocate(bscoef_density(n1,n2,n3))
    allocate(xknot(nxknot))
    allocate(yknot(nyknot))
    allocate(zknot(nzknot))
    call dbs3ines(xx,yy,zz,V,n1,n2,n3,nxknot,nyknot,nzknot,bscoef_V,xknot,yknot,zknot)
    call dbs3ines(xx,yy,zz,density,n1,n2,n3,nxknot,nyknot,nzknot,bscoef_density,xknot,yknot,zknot)
    if (DFT) then
       call grad_nabla_interp(n1,n2,n3,nxknot,nyknot,nzknot,xknot,yknot,zknot,xx,yy,zz,nAtom,x,y,z)
    end if
    write(6,*)'done  interpolations'
    !-----------------
    !*COMM* find charge to subtract form singularity at nucleus
    !do i=1,nAtom
    !val=dbs3vl(x(i)-epsilon,y(i)-epsilon,z(i)-epsilon,order,order,order,xknot,yknot,zknot,n1,n2,n3,bscoef)
    !charge_nucl(i)=val*(abs(position(i)- Sqrt((x(i)-epsilon)**2+(y(i)-epsilon)**2+(z(i)-epsilon)**2)))
    !write(6,*)charge_nucl(i)
    !end do
    !-----------------
    !write(6,*)'alpha='
    !read(5,*)alpha
    !XXXXXX
    !stop
    !XXXXXX
    deallocate(V,STAT=test)
    deallocate(density)
    !deallocate(xx)
    !deallocate(yy)
    !deallocate(zz)
    ! New deallocations
    deallocate(r_glo_avg)
    deallocate(theta_glo_avg)
    deallocate(phi_glo_avg)
    !write(6,*)'alpha=',test
    !read(5,*)alpha
deallocate(Atomic_num,Charge)
!deallocate(x,y,z)
deallocate(theta_d,raggio_d,rho_d,z_d)
    write(6,*)'end V_conversion'
  end subroutine V_conversion



  subroutine grad_nabla_interp(n1,n2,n3,nxknot,nyknot,nzknot,xknot,yknot,zknot,xx,yy,zz,nAtom,x,y,z)
    use bspline , only : dbs3in,dbs3vl

    implicit none
    integer(kind=i4b), intent(in) :: n1,n2,n3,nxknot,nyknot,nzknot,nAtom
    real(kind=dbl), intent(in) :: xx(n1),yy(n2),zz(n3),xknot(nxknot),yknot(nyknot),zknot(nzknot),x(nAtom),y(nAtom),z(nAtom)
    real(kind=dbl), dimension(:,:,:),allocatable :: grad,nabla
    integer (kind=i4b) :: nx,ny,i1,i2,i3,i,test,na
    real(kind=dbl) :: val_grad,val_nabla,xt,yt,zt
    !-------------------------
    !1st part
    nx=13
    ny=14
    open(unit=nx,file='grad.dat',status='old',action='read')
    open(unit=ny,file='nabla.dat',status='old',action='read')
    do i=1,6+nAtom
       read(ny,*)
       read(nx,*)
    end do
    allocate(grad(n1,n2,n3),STAT=test)
    write(6,*)test
    allocate(nabla(n1,n2,n3),STAT=test)
    write(6,*)test
    do i1=1,n1
       do i2=1,n2
          read(nx,900)( grad(i1,i2,i3),i3=1,n3)
          !write(3009,900)(grad(i1,i2,i3),i3=1,n3)
          read(ny,900)( nabla(i1,i2,i3),i3=1,n3)
       end do
    end do
    do i1=1,n1
       do i2=1,n2
          do i3=1,n3
             !grad(i1,i2,i3)=sqrt(grad(i1,i2,i3))
             !if (abs(yy(i2)).le..1d0) then
             !   if (abs(xx(i1)).le..05d0) then
             !      write(3007,*)grad(i1,i2,i3),zz(i3)
             !      write(3008,*)nabla(i1,i2,i3),zz(i3)
             !   end if
             !end if
             do na=1,nAtom
                nabla(i1,i2,i3) = nabla(i1,i2,i3) * ((xx(i1)-x(na))**2+(yy(i2)-y(na))**2+(zz(i3)-z(na))**2)
             end do
             do na=1,nAtom
                grad(i1,i2,i3) = grad(i1,i2,i3) * ((xx(i1)-x(na))**2+(yy(i2)-y(na))**2+(zz(i3)-z(na))**2)
             end do
             !if (abs(yy(i2)).le..1d0) then
             !   if (abs(xx(i1)).le..05d0) then
             !      write(3009,*)grad(i1,i2,i3),zz(i3)
             !      write(3010,*)nabla(i1,i2,i3),zz(i3)
             !   end if
             !end if
          end do
       end do
    end do
    allocate(bscoef_grad(n1,n2,n3))
    allocate(bscoef_nabla(n1,n2,n3))
    call dbs3ines(xx,yy,zz,grad,n1,n2,n3,nxknot,nyknot,nzknot,bscoef_grad,xknot,yknot,zknot)
    call dbs3ines(xx,yy,zz,nabla,n1,n2,n3,nxknot,nyknot,nzknot,bscoef_nabla,xknot,yknot,zknot)
    !-----------------------------------------
    ! Verify potential
    !do i1=1,300
    !   do i2=1,300
    !      yt=0.d0
    !      xt=-5.d0+i1*.03d0
    !      zt=-5.d0+i2*.03d0
    !      val_grad=dbs3vl(xt,yt,zt,5,5,5,xknot,yknot,zknot,n1,n2,n3,bscoef_grad)
    !      val_nabla=dbs3vl(xt,yt,zt,5,5,5,xknot,yknot,zknot,n1,n2,n3,bscoef_nabla)
    !      do na=1,nAtom
    !         val_grad=val_grad / ((xt-x(na))**2+(yt-y(na))**2+(zt-z(na))**2)
    !      end do
    !      do na=1,nAtom
    !         val_nabla=val_nabla / ((xt-x(na))**2+(yt-y(na))**2+(zt-z(na))**2)
    !      end do

    !      !Check interpolation
    !      write(30,*)xt,zt,val_grad
    !      write(31,*)xt,zt,val_nabla
    !      !End check
    !   end do
    !end do
    !-------------------------------------------
    !write(109,*)bscoef_grad
    deallocate(grad)
    deallocate(nabla)
9000 format(1A40)
900 format(6E13.5)
    write(6,*)'end grad_nabla_interp'
  end subroutine grad_nabla_interp
  !-------------------------

  subroutine grad_nabla_integr_prep(n,max_num,sing_sectors,cube_num_pts_sing)
    use   Calc_func,only :  interp_grad,interp_grad_sing,interp_nabla,interp_nabla_sing
    implicit none
    integer(kind=i4b), intent(in) :: n,max_num,sing_sectors,cube_num_pts_sing
    integer(kind=i4b) :: i,j
    write(6,*)'grad_nabla_integr_prep'
    allocate(interp_grad(n,max_num))
    allocate(interp_nabla(n,max_num))
    do i=1,n
       do j=1,max_num
          interp_grad(i,j)=0.0d0
          interp_nabla(i,j)=0.0d0
       end do
    end do
    allocate(interp_grad_sing(sing_sectors,cube_num_pts_sing))
    allocate(interp_nabla_sing(sing_sectors,cube_num_pts_sing))
    do i=1,sing_sectors
       do j=1,cube_num_pts_sing
          interp_grad_sing(i,j) = 0.0d0
          interp_nabla_sing(i,j) = 0.0d0
       end do
    end do
  end subroutine grad_nabla_integr_prep

  subroutine integral_release_memory(max_num,sing_sectors,cube_num_pts_sing)
    use   Calc_func
    !,only :  interp_grad,interp_grad_sing,interp_nabla,interp_nabla_sing
    implicit none
    integer(kind=i4b), intent(in) :: max_num,sing_sectors,cube_num_pts_sing
    write(6,*)'grad_nabla_integr_prep release memory'
    if (DFT) then
       deallocate(interp_grad)
       deallocate(interp_nabla)
       deallocate(interp_grad_sing)
       deallocate(interp_nabla_sing)
    end if
    deallocate(interp_V)
    deallocate(interp_density)
    deallocate(vect_sing)
    deallocate(interp_V_sing)
    deallocate(interp_density_sing)
    deallocate(val_array)
    deallocate(val_array_sing) 
    deallocate(bas_tab)
    deallocate(bas_tab_sing)
    deallocate(der_tab)
    deallocate(der_tab_sing)
    deallocate(mem_tab)

  end subroutine integral_release_memory


end module V_setup
