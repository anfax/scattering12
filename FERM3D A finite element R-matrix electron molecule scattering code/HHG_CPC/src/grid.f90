!!!*************************************************************
! 文件/File: grid.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: grid.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


subroutine grid
  use Brep
  use control
  use calc_func, only : polar,polar_2,r_cutoff,Energy_last_bound
  use nrtype , only : Pi,dbl,i4b
  !  use bsp_info,only : bsppid, bspnprocs
  use Open_information, only : max_l 
  use function_integr_tab,only : r_glo_in_tab,r_glo_in_t,theta_glo_in_tab,theta_glo_in_t
  use Solve, only : charge=> Z
  implicit none 
  real(kind=dbl) :: r,phi,phitemp,rtemp,theta,thetatemp,x_coord,y_coord,z_coord
  integer(kind=i4b) :: zone,nr,nphi,globnode,element,i,j,ntheta
  integer(kind=i4b) :: dim
  real(kind=dbl),dimension(:),allocatable ::Delta_rp,Delta_phip,Delta_thetap,delta_r,delta_theta,delta_phi
  integer(kind=i4b),dimension(:),allocatable :: nr_ztemp,nphi_ztemp,ntheta_ztemp
  real(kind=dbl),dimension(:),allocatable :: nr_zp,nphi_zp,ntheta_zp

  !-----------------------------------------------------------------------
  ! This subroutine makes the grids, calculates the vertices (r,phi,theta) and 
  ! makes the elements, then makes the global indexing of the functions.
  ! Variables:
  ! zone = identifies the parts of the domain in which there are different
  !        grids
  ! nr,nphi = indices of number of sectors in r,phi for each zone
  ! nr_z,nphi_z = number of sectors in r,phi for each zone
  ! globnode,element = indices for the physical nodes and the sectors
  ! i,j = indices
  ! dim = max dimension of the arrays of the grid arrays
  ! Delta_rp,Delta_phip,nr_ztemp,nphi_ztemp,nr_zp,nr_phip = temporary Delta_r,
  !                            Delta_phi,nr_ztemp,nphi_ztemp,nr_zp,nr_phip
  !-----------------------------------------------------------------------
  ! initialization
  !R0=box size
  !*COMM******************
  open(unit=10,file='input_grid.dat',status='old',action='read')
  open(unit=12,file='input_control.dat',status='old',action='read')
  read(12,*)output_file
  read(12,*)calculation_type,option_wf
  read(12,*)exchange,max_l
  read(12,*)molecule,charge

  read(12,*)polar
  read(12,*)polar_2
  read(12,*)r_cutoff
  read(12,*)Energy_last_bound
  read(12,*)shift
  close(12)
  !*COMM******************
  ! Write to std output the input flags
  write(6,*)'-----------------------------------------'
  write(6,*)'CALCULATION FLAGS'
  write(6,*)'output file= ',output_file
  write(6,*)'calculation type= ',calculation_type,'wf= ',option_wf
  write(6,*)'exchange= ',exchange
  write(6,*)'maximum continuum angular momentum', max_l
  write(6,*)'molecule type= ',molecule,'charge= ',charge
  write(6,*)'polarizability alpha_0= ',polar
  write(6,*)'polarizability alpha_2= ',polar_2
  write(6,*)'r_cutoff= ',r_cutoff
  write(6,*)'Energy last bound state= ',Energy_last_bound
  write(6,*)'shift in z-direction= ',shift
  write(6,*)'-----------------------------------------'

  !*COMM******************
  read(10,*)R0
  write(6,*)'R0=',R0
  ! z=number of zones of the grid with different meshes 
  read(10,*)z

  ! Input checking
write(6,*)calculation_type,'mostro'
  if ((calculation_type.ne.'scatter').and.(calculation_type.ne.'bound')) then 
    write(6,*)'error, the calculation type is wrong, needs to be <scatter> or <bound>'
    stop
  end if
  if ((option_wf.ne.'wf').and.(option_wf.ne.'nowf')) then
    write(6,*)'error, the wavefunction printing option is wrong, needs to be <wf> or <nowf>'
    stop
  end if
  if (exchange.ne.'slater') then
    write(6,*)'error, the exchange variable has to be <slater>'
    stop
  end if
  if (max_l.lt.1) then
    write(6,*)'error, the maximum angular momentum has to be max_l>0'
    stop
  end if
  if ((molecule.ne.'neutral').and.(molecule.ne.'ion')) then
    write(6,*)'error, the molecule type is wrong, needs to be <neutral> or <ion>'
    stop
  end if
  if (charge.lt.0) then
    write(6,*)'error, the molecular charge needs to be greater or equal to zero'
    stop
  end if
  if ((polar.lt.0.d0).or.(polar_2.lt.0.d0)) then
    write(6,*)'error, the polarizability coefficients need to be greater or equal to zero'
    stop
  end if
  if ((r_cutoff.le.0.d0).or.(r_cutoff.ge.R0)) then
    write(6,*)'error, the cutoff radius needs to be greater than zero and lower than R0, the R-matrix radius'
    stop
  end if
  if (Energy_last_bound.lt.0) then
    write(6,*)'error, the ionization energy needs to be greater or equal to zero'
    stop
  end if
  if (R0.le.0.d0) then
    write(6,*)'error, the R-matrix radius needs to be greater than zero'
    stop
  else if (R0.le.5.d0) then
    write(6,*)'warning, the R-matrix radius is',R0,'a very small value'
  end if

  !----------------------------------------------------------------------------

9000 format(2a10)
9001 format(1a10,1f13.5)
  !dimensioning of node and element arrays
  dim=100000
  !-----------------------------------------------------------------------
  !1st allocations
  allocate(r_in(z))
  allocate(r_out(z))
  allocate(phi_in(z))
  allocate(phi_out(z))
  allocate(theta_in(z))
  allocate(theta_out(z))

  allocate(Delta_r(z))
  allocate(Delta_phi(z))
  allocate(Delta_theta(z))
  allocate(nr_z(z))
  allocate(nphi_z(z))
  allocate(ntheta_z(z))

  allocate(Delta_rp(z))
  allocate(Delta_phip(z))
  allocate(Delta_thetap(z))
  allocate(nr_zp(z))
  allocate(nphi_zp(z))
  allocate(ntheta_zp(z))

  !-----------------------------------------------------------------------


  !parameters of the zones
  !---------------------------
  do zone=1,z

     read(10,*)r_in(zone),r_out(zone),phi_in(zone),phi_out(zone),theta_in(zone),theta_out(zone),Delta_rp(zone),Delta_phip(zone),Delta_thetap(zone) 
9090 format(9f13.8)
     nr_zp(zone)=(r_out(zone)-r_in(zone))/Delta_rp(zone)
     nphi_zp(zone)=(phi_out(zone)-phi_in(zone))/Delta_phip(zone)
     ntheta_zp(zone)=(theta_out(zone)-theta_in(zone))/Delta_thetap(zone)    
     !$     write(29,*)nphi_zp(zone),nr_zp(zone)
     nr_z(zone)=anint(nr_zp(zone))
 if (nr_z(zone).le.0) nr_z(zone)=1

     nphi_z(zone)=anint(nphi_zp(zone))
 if (nphi_z(zone).le.0) nphi_z(zone)=1

     ntheta_z(zone)=anint(ntheta_zp(zone))
 if (ntheta_z(zone).le.0) ntheta_z(zone)=1

     Delta_r(zone)=(r_out(zone)-r_in(zone))/nr_z(zone)
     Delta_phi(zone)=(phi_out(zone)-phi_in(zone))/nphi_z(zone)
     Delta_theta(zone)=(theta_out(zone)-theta_in(zone))/ntheta_z(zone)        
     Delta_phi(zone)=Delta_phi(zone)*Pi
     phi_in(zone)=phi_in(zone)*Pi  
     phi_out(zone)=phi_out(zone)*Pi
     Delta_theta(zone)=Delta_theta(zone)*Pi
     theta_in(zone)=theta_in(zone)*Pi
     theta_out(zone)=theta_out(zone)*Pi 
     phi_sectors=2*Pi/Delta_phi(zone)
!!$write(29,*)Delta_r(zone),Delta_phi(zone)
  end do


  !-----------------------------------------------------------------------
  ! Grid checking --> See if more needed
  do zone=1,z
    if ((r_out(zone).lt.R0).and.(theta_out(zone).lt.Pi).and.(phi_out(zone).lt.2.d0*Pi)) then
     do nr=1,z
       if ((r_out(zone).eq.r_in(nr)).and.&
        & (theta_out(zone).eq.theta_in(nr)).and.&
        & (phi_out(zone).eq.phi_in(nr))) then
           exit
       else if (nr.eq.z) then
           write(6,*)'grid error, check the input grid'
           stop
       end if
     end do
    else if ((r_out(zone).lt.R0).and.(theta_out(zone).lt.Pi).and.(phi_out(zone).eq.2.d0*Pi)) then
     do nr=1,z
       if ((r_out(zone).eq.r_in(nr)).and.&
        & (theta_out(zone).eq.theta_in(nr)).and.&
        & (phi_in(nr).eq.0.d0)) then
           exit
       else if (nr.eq.z) then
           write(6,*)'grid error, check the input grid'
           stop
       end if
     end do
    end if
  end do
  if (r_out(z).lt.R0) then
      write(6,*)'grid error, check the input grid'
      stop
  end if
  !-----------------------------------------------------------------------

  deallocate(nr_zp)
  deallocate(nphi_zp)
  deallocate(ntheta_zp)
  deallocate(Delta_rp)
  deallocate(Delta_phip)
  deallocate(Delta_thetap)
  !-----------------------------------------------------------------------
90 format(6d14.7)
  globnode=0


  allocate(nr_ztemp(z))
  allocate(nphi_ztemp(z))
  allocate(ntheta_ztemp(z))
  allocate(r_glo(dim))
  allocate(phi_glo(dim))
  allocate(theta_glo(dim))
  allocate(zone_flag(dim))
  allocate(boundary(dim))
  !----------------------------------------------------------------------- 
  r=0.0d0
  phi=0.0d0
  !-----------------------------------------------------------------------

  !nodes making

  !-----------------------------------------------------------------------
  !$  write(17,*)R0
  do zone=1,z
     r=r_in(zone)
     if (r_out(zone).ge.(R0-epsilon_custom)) then
        nr_z(zone)=nr_z(zone)+1
     end if
     !-----------------------------------------------------------------------
     !Needed to count nodes at 2*Pi as different for rectangular grid
     !  if((phi_out(zone).eq.2)) then
     !    nphi_z(zone)=nphi_z(zone)+1
     !  end if
     !-----------------------------------------------------------------------
     if((theta_out(zone).ge.(Pi-epsilon_custom))) then
        ntheta_z(zone)=ntheta_z(zone)+1
     end if

     do nr=1,nr_z(zone)
        phi=phi_in(zone)
        do nphi=1,nphi_z(zone)
           theta=theta_in(zone)
           do ntheta=1,ntheta_z(zone)
              globnode=globnode+1
              r_glo(globnode)=r
              phi_glo(globnode)=phi
              theta_glo(globnode)=theta
              !              write(102,*)theta_glo(globnode),globnode
              !      write(13,*)r_glo(globnode),R0
              if (r.ge.(R0-epsilon_custom)) then
                 boundary(globnode)='up'

              else if (r_glo(globnode).eq.0) then
                 boundary(globnode)='down'

              else if (phi_glo(globnode).eq.0) then
                 boundary(globnode)='left'
              else if (phi_glo(globnode).ge.(2.0d0*Pi-epsilon_custom)) then
                 !      else if (phi_glo(globnode).ge.Pi*2.0d0-(1.0d-4)) then
                 boundary(globnode)='right'
                 !$                 write(101,*)'giorgia'
              else if (theta_glo(globnode).ge.(Pi-epsilon_custom)) then
                 boundary(globnode)='rear'
              else if  (theta_glo(globnode).eq.0) then
                 boundary(globnode)='front'
              else
                 boundary(globnode)='bulk'
              end if

              theta=theta+Delta_theta(zone)
              !write(8,*)r_glo(globnode),phi_glo(globnode),theta_glo(globnode),globnode,boundary(globnode)
              !write(9,*)r_glo(globnode),phi_glo(globnode),theta_glo(globnode)
           end do
           phi=phi+Delta_phi(zone)
           !                 write(7,*)r_glo(globnode),phi_glo(globnode),globnode,boundary(globnode)
        end do
        r=r+Delta_r(zone)
     end do
  end do
  !XX!
  !------------------------------------
  !elements making, there are 8 nodes for each physical node localized in space,
  !to account for the 4 basis functions for each dimension   

  allocate(nodes(dim,8))
  allocate(PatchVertex(dim,64))
  r=0.0d0
  phi=0.0d0
  theta=0.0d0
  do i=1,dim
     do j=1,8
        nodes(i,j)=0
     end do
     do j=1,64
        PatchVertex(i,j)=0
     end do
  end do
  element=0
  p=globnode
  allocate(delel_r(p))
  allocate(delel_theta(p))
  allocate(delel_phi(p))
  do zone=1,z
     if (r_out(zone).ge.(R0-epsilon_custom)) then
        nr_z(zone)=nr_z(zone)-1
     end if
     !---------------------------------------------
     ! Continuity 0-2*Pi on the Phi variable
     !---------------------------------------------
     !---------------------------------------------
     if((phi_out(zone).ge.(2*Pi-epsilon_custom))) then
        !    nphi_z(zone)=nphi_z(zone)-1
     end if
     if((theta_out(zone).ge.(Pi-epsilon_custom))) then
        ntheta_z(zone)=ntheta_z(zone)-1
     end if
     r=r_in(zone)


     !---------------------------------------------
     !---------------------------------------------
     do nr=1,nr_z(zone)
        phi=phi_in(zone)
        do nphi=1,nphi_z(zone)
           theta=theta_in(zone)
           do ntheta=1,ntheta_z(zone)
              rtemp=r
              phitemp=phi
              thetatemp=theta
              element=element+1
              zone_flag(element)=zone
              delel_r(element)=Delta_r(zone)
              delel_theta(element)=delta_theta(zone)
              delel_phi(element)=delta_phi(zone)
              do i=1,8
                 phi=phitemp
                 r=rtemp
                 theta=thetatemp
                 if (i.eq.2) then
                    theta=theta+Delta_theta(zone)
                 else if (i.eq.3) then
                    phi=phi+Delta_phi(zone)
                 else if (i.eq.4) then
                    theta=theta+Delta_theta(zone)
                    phi=phi+Delta_phi(zone)
                 else if (i.eq.5) then
                    r=r+Delta_r(zone)
                 else if (i.eq.6) then
                    r=r+Delta_r(zone)
                    theta=theta+Delta_theta(zone)
                 else if (i.eq.7) then
                    r=r+Delta_r(zone)
                    phi=phi+Delta_phi(zone)
                 else if (i.eq.8) then
                    theta=theta+Delta_theta(zone)
                    phi=phi+Delta_phi(zone)
                    r=r+Delta_r(zone)
                 end if
                 globnode=0


                 do globnode=1,p

                    !--------------------------------------------------------------------------------
                    ! Imposing continuity at each element boundary by numbering of the physical nodes
                    ! PatchVertex array makes the connection between the 64 local nodes in each
                    ! element (used to do integrals)
                    ! and the global nodes that will be used to construct the matrices
                    !--------------------------------------------------------------------------------

                    if (phi.ge.(Pi*2.0d0-epsilon_custom)) then
                       if(((r_glo(globnode).ge.(r-epsilon_custom)).and.(r_glo(globnode).le.(r+epsilon_custom)))&
                            &                 .and.(phi_glo(globnode).le.epsilon_custom).and.((theta_glo(globnode).ge.(theta-epsilon_custom))&
                            &                .and.(theta_glo(globnode).le.(theta+epsilon_custom)))) then
                          nodes(element,i)=globnode
                          do j=1,8
                             PatchVertex(element,8*i-(8-j))=8*globnode-(8-j)
                          end do
                       end if

                    else if(((r_glo(globnode).ge.(r-epsilon_custom)).and.(r_glo(globnode).le.(r+epsilon_custom)))&
                         &                     .and.((phi_glo(globnode).ge.(phi-epsilon_custom)).and.(phi_glo(globnode).le.(phi+epsilon_custom)))&
                         &                     .and.((theta_glo(globnode).ge.(theta-epsilon_custom)).and.(theta_glo(globnode).le.(theta+epsilon_custom))) ) then
                       nodes(element,i)=globnode
                       do j=1,8
                          PatchVertex(element,8*i-(8-j))=8*globnode-(8-j)
                       end do

                    end if
                 end do
              end do
              theta=thetatemp
              phi=phitemp
              r=rtemp
              theta=theta+Delta_theta(zone)
           end do
           phi=phi+Delta_phi(zone)
        end do
        r=r+Delta_r(zone) 
     end do
  end do
  n=element
  p_Patch=8*p
  write(6,*) 'done grid, elements, Patches=',n,p

  !---------------
  !*COMM* Element checking
  !        do element=1,n
  !          do i=1,8
  !            write(451,*)r_glo(nodes(element,i)),element,i
  !          end do
  !        end do
  !----------------

  !----------------
  !*COMM* coordinates for wavefunction checking
  allocate(coord_x(p))
  allocate(coord_y(p))
  allocate(coord_z(p))
  coord_x=0.d0
  coord_y=0.d0
  coord_z=0.d0
plotpoints=0
allocate(nodes_closed_index(p))
  do globnode=1,p
!     if (boundary(globnode).ne.'up') then
   if (r_glo(globnode).lt.r_in(z)) then
        plotpoints=plotpoints+1
        coord_x(globnode)=r_glo(globnode)*sin(theta_glo(globnode))*cos(phi_glo(globnode))
        coord_y(globnode)=r_glo(globnode)*sin(theta_glo(globnode))*sin(phi_glo(globnode))
        coord_z(globnode)=r_glo(globnode)*cos(theta_glo(globnode))
        write(3888,*)coord_x(globnode),coord_y(globnode),coord_z(globnode)
        nodes_closed_index(globnode)=plotpoints
     end if
  end do
  !        return
  !----------------
  deallocate(nr_ztemp)
  deallocate(nphi_ztemp)
  deallocate(ntheta_ztemp)
  deallocate(zone_flag)
  !XXXXXXX
  call r_glo_in_t
  call theta_glo_in_t
deallocate(r_in,r_out,Delta_r,Delta_phi,Delta_theta)
deallocate(nr_z,nphi_z,ntheta_z)
deallocate(theta_in,theta_out,phi_in,phi_out)
!deallocate(phi_glo,r_glo,theta_glo)
  !XXXXXXX
end subroutine grid


real(kind=8) function r_glo_in(element)
  use Brep
  use nrtype, only : dbl,i4b
  implicit none
  integer(kind=i4b) :: l
  integer(kind=i4b),intent(in) :: element
  integer(kind=i4b),parameter :: first=1
  real(kind=dbl) :: r_glo_in_temp
  ! This subroutine calculates the initial radius value for each element/sector
  r_glo_in_temp= r_glo(nodes(element,first))
  do l=1,8
     if (r_glo(nodes(element,l)).lt.r_glo_in_temp) then
        r_glo_in_temp=r_glo(nodes(element,l)) 
     end if
  end do
  r_glo_in=r_glo_in_temp
  return
end function r_glo_in


subroutine r_glo_average
  use Brep
  use nrtype, only : i4b
  implicit none
  integer(kind=i4b) :: element
  integer(kind=i4b),parameter :: first=1,last=8
  !This subroutine calculates the value of the radius at the center of the sector

  allocate(r_glo_avg(n))
  !  write(381,*)r_glo
  do element=1,n
     r_glo_avg(element)=(r_glo(nodes(element,first))+r_glo(nodes(element,last)))/2.d0
     !     write(378,*)r_glo(nodes(element,1))
  end do
end subroutine r_glo_average


subroutine theta_glo_average
  use Brep
  use nrtype, only : i4b
  implicit none
  integer(kind=i4b) :: element
  integer(kind=i4b),parameter :: first=1,last=8
  !This subroutine calculates the value of  theta at the center of the sector

  allocate(theta_glo_avg(n))
  !  write(381,*)r_glo
  do element=1,n
     theta_glo_avg(element)=(theta_glo(nodes(element,first))+theta_glo(nodes(element,last)))/2.d0
     !     write(378,*)r_glo(nodes(element,1))
  end do
end subroutine theta_glo_average

real(kind=8) function theta_glo_in(element)
  use nrtype, only : dbl,i4b
  use Brep
  implicit none
  integer(kind=i4b) :: l 
  integer(kind=i4b),intent(in):: element
  integer(kind=i4b),parameter :: first=1
  real(kind=dbl) :: theta_glo_in_temp
  ! This subroutine calculates the initial theta value for each element/sector
  theta_glo_in_temp= theta_glo(nodes(element,first))
  do l=1,8
     if (theta_glo(nodes(element,l)).lt.theta_glo_in_temp) then
        theta_glo_in_temp=theta_glo(nodes(element,l))
     end if
  end do
  theta_glo_in=theta_glo_in_temp
  return
end function theta_glo_in




real(kind=8) function phi_glo_in(element)
  use nrtype, only :dbl,i4b
  use Brep,only : phi_glo,nodes,epsilon_custom
  implicit none
  integer(kind=i4b) :: l 
  integer(kind=i4b),intent(in):: element
  integer(kind=i4b),parameter :: first=1
  real(kind=dbl) :: phi_glo_in_temp
  !write(332,*)element
  phi_glo_in_temp= phi_glo(nodes(element,first))
  do l=1,8
     if (phi_glo(nodes(element,l)).ge.epsilon_custom) then
        if (phi_glo(nodes(element,l)).lt.phi_glo_in_temp) then
           phi_glo_in_temp=phi_glo(nodes(element,l)) 
        end if
     end if
  end do
  phi_glo_in=phi_glo_in_temp
end function phi_glo_in


subroutine phi_glo_average
  use Brep
  implicit none
  integer :: element
  integer,parameter :: first=1,last=8
  !This subroutine calculates the value of the phi angle at the center of the sector
  allocate(phi_glo_avg(n))
  do element=1,n
     phi_glo_avg(element)=(phi_glo(nodes(element,first))+phi_glo(nodes(element,last)))/2.d0
  end do
end subroutine phi_glo_average


subroutine potential
  use Brep, only : r_glo_avg, theta_glo_avg,phi_glo_avg,U,n 
  use nrtype, only : Pi,dbl,i4b
  use V_setup, only : V_conversion
  implicit none
  real(kind=dbl),parameter :: k=1.0d0
  real(kind=dbl),dimension(:),allocatable :: charge_density
  integer(kind=i4b) :: element,num_gauss,i
  integer(kind=i4b),parameter :: nx=11
  !  allocate(U(n))
  !  allocate(charge_density(n))
  ! This subroutine makes the calls for the reading and interpolation of the potential

  call r_glo_average
  call theta_glo_average
  call phi_glo_average
  !						
  !	unit name: nx=11 for the potential	
  !						
  call V_conversion(n,nx)
  !call charge(charge_density)
end subroutine potential
