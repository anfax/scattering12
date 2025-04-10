!!!*************************************************************
! 文件/File: prog_CM.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: prog_CM.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


program CM

  ! Program to do gridding for FERM3D.x, on the basis of the position of the
  ! atoms in the molecule.
  ! Some fixed points are: R0= R-matrix radius, R_cutoff+0.5d0, theta=Pi, phi=2Pi
  ! There are  an internal and external region with different radial grids

  integer :: i,j,f,k,natom,n,signo,kk,counter,kkk
  integer,dimension(3) ::f_store 
  real*8,dimension(:),allocatable :: x,y,z,charge,delta,tempa,tempb,r
  real*8,dimension(:,:),allocatable :: cross_product,coord,coord_a
  real*8 :: XCM,YCM,ZCM,M_tot,dummy,temp1,temp2,shiftx,shifty,shiftz,R_cutoff
  real*8,parameter,dimension(3) :: thresh=(/1.d-5,1.d-5,3.d-5/)
  real*8,parameter ::  Pi=3.141592653589793238462643383279502884197d0,&
       &ang_to_au=.5291d0,&
       &epsilon_custom=1.d-5
  real*8 :: R0=12d0,delta_theta_max=0.0825d0,delta_phi_max=0.0625d0,&
       &delta_r_max=0.125d0,delta_external_radius=.5d0
  logical :: stat,rotation=.TRUE.,control_zero_phi=.FALSE.
  logical, parameter :: angstrom=.FALSE.,CentMass=.TRUE.
  integer, parameter :: sumdim=30
  ! Coordinate 1=r, 2=phi, 3=theta
  ! Read CM.dat information, grid data
  open(unit=7,file='CM.dat')
  read(7,*)R0
  read(7,*)delta_theta_max
  read(7,*)delta_phi_max
  read(7,*)delta_r_max
  read(7,*)delta_external_radius
  close(7)
  ! Read input_control.dat, input information
  open(unit=7,file='input_control.dat')
  do i=1,6
     read(7,*)
  end do
  read(7,*)R_cutoff
  read(7,*)
  read(7,*)shiftz
  do i=1,3
     read(7,*)
  end do
  read(7,*)shifty
  read(7,*)rotation
  read(7,*)
  read(7,*)shiftx
  close(7)
  write(6,*)'input data',R_cutoff,shiftz,shiftx,shifty

  ! Read molecular data from pot.dat
  open(unit=8,file='pot.dat')
  do i=1,2
     read(8,*)
  end do
  read(8,*)natom
  write(6,*)natom
  do i=1,3
     read(8,*)
  end do
  allocate(x(natom),y(natom),z(natom),charge(natom))
  allocate(coord(0:natom+sumdim,3),coord_a(0:natom+sumdim,3))
  allocate(delta(3))

  do i=1,natom
     read(8,*)temp1,charge(i),x(i),y(i),z(i)
  end do
  close(8)
  M_tot=0.d0
  XCM=0.d0;YCM=0.d0;ZCM=0.d0
  do i=1,natom
     M_tot=M_tot+charge(i)
     XCM=XCM+x(i)*charge(i)
     YCM=YCM+y(i)*charge(i)
     ZCM=ZCM+z(i)*charge(i)
  end do

  XCM=XCM/M_tot
  YCM=YCM/M_tot
  ZCM=ZCM/M_tot

  write(6,*)XCM,YCM,ZCM
  !--------------------
  ! Shift according to input_control.dat
  XCM=shiftx
  YCM=shifty
  ZCM=shiftz
  !--------------------
  ! Center of mass translation

  do i=1,natom
     if (angstrom.and.CentMass) then
        x(i)=(x(i)-XCM)/ang_to_au
        y(i)=(y(i)-YCM)/ang_to_au
        z(i)=(z(i)-ZCM)/ang_to_au
     else if (CentMass) then
        x(i)=(x(i)-XCM)
        y(i)=(y(i)-YCM)
        z(i)=(z(i)-ZCM)
     else
        x(i)=(x(i)-XCM)
        y(i)=(y(i)-YCM)
        z(i)=(z(i)-ZCM)
        write(6,*)angstrom,CentMass
     end if
     !write(9,*)x(i),y(i),z(i)
  end do
  do i=1,natom
     !---------------------
     ! Pseudo rotation of the axes
     if (rotation) then
        temp2=y(i)
        y(i)=z(i)
        z(i)=temp2
     end if
     !----------------------
  end do
  ! Polar coordinates transformation
  ! There will be a sector boundary passing through each atomic nucleus; 
  !coord(i,1)=r coord for atom i, coord(i,2)=phi coord for atom i, coord(i,3)=theta coord for atom i

  do i=1,natom
     coord(i,1)=sqrt(x(i)**2+y(i)**2+z(i)**2)
     if (coord(i,1).gt. 0.0d0+epsilon_custom) then
        coord(i,3)=acos(z(i)/coord(i,1))
        !        if ((abs(x(i)).gt.1.d-12)) then
        if ((coord(i,3).ge.epsilon_custom).and.(coord(i,3).le.Pi-epsilon_custom)) then
           if (abs (x(i)/(coord(i,1)*sin(coord(i,3)))).le.1.d0) then
              coord(i,2)=acos(x(i)/(coord(i,1)*sin(coord(i,3))))
           else if ( (x(i)/(coord(i,1)*sin(coord(i,3)))).lt.-1.d0) then
              coord(i,2)=Pi
           else if ( (x(i)/(coord(i,1)*sin(coord(i,3)))).gt.1.d0) then
              coord(i,2)=0.d0
           end if
           !     if (abs(x(i)).gt.1.d-12) then
           !     phi(i)=atan(y(i)/x(i))
           if (y(i).lt.0.d0) coord(i,2)=2.d0*Pi-coord(i,2)
        else
           if (y(i).gt.0.d0) then
              coord(i,2)=Pi/2.d0
           else if (y(i).lt.0.d0) then
              coord(i,2)=3*Pi/2.d0 
           else if (y(i).eq.0.d0) then
              coord(i,2)=0.d0
           end if
        end if
     else
        coord(i,3)=0.d0
        coord(i,2)=0.d0
     end if
     write(10,*)coord(i,1),coord(i,3),coord(i,2)
  end do

  !Now for each atom, need a zone boundary at k-shell radius for r, theta, phi: 
  !one at nucleus plus k-shell radius, one at nucleus minus k-shell radius.
  !tempa stores coords
  allocate(tempa(0:2*natom+sumdim),tempb(0:2*natom+sumdim),r(0:2*natom+sumdim))
  ! Choose points
  !put coords of nuclei into grid, as well as other points which we know have 
  !to be in the grid: zero(for all points)
  !in r: 0, r_cutoff, r_0 (r_cutoff is where you switch between short-range and long-range treatment of polarization)
  !in theta: 0, Pi
  !in phi: 0, 2*Pi
  do kk=1,3!k ranges over dimensions
  coord_a(:,kk)=0.d0
     counter=2!counter counts how many points have already been put into each dimension
     if (kk.eq.1) then!r dimension
        coord_a(counter,1)=R_cutoff+0.5d0!+1.5d0 might be a better value, because
        counter=counter+1
        coord_a(counter,1)=R0
     else if (kk.eq.3) then!theta dimension
        coord_a(counter,3)=Pi
     else if (kk.eq.2) then!phi dimension
        coord_a(counter,2)=2.d0*Pi
     end if
     stat=.TRUE.
     f=counter
     !now put coords of nuclei into array, provided they're not already in the array via earlier process
     do i=1,natom
        do j=1,f
           stat=.FALSE.
           if (abs((coord(i,kk)-coord_a(j,kk))).gt.thresh(kk)) then
              stat=.TRUE.
           else
              stat=.FALSE.
              exit
           end if
           if (j==f.and.STAT) then
              coord_a(f+1,kk)=coord(i,kk)
           end if
        end do
        if (stat) then
           f=f+1
           stat=.FALSE.
        end if
     end do
     do i=0,f
        write(11,*)coord_a(i,kk),i,kk
     end do
     f_store(kk)=f
!+counter
  end do

  !now put points at k-shell radii of each atom in each dimension
  !-------------------------------------------------
  ! K-shell radius
  do kk=1,3!kk ranges over dimensions
     f=f_store(kk)
     tempb(:)=0.d0;tempa(:)=0.d0
     tempb(0:f_store(kk))=coord_a(0:f_store(kk),kk)
     tempa(:)=coord(:,kk)
write(35,*)'aa',tempa(:)
write(35,*)'bb',tempb(:)
     if (kk.eq.1) then!r dimension
        r(:)=1.d0
        kkk=3
     else 
        r(:)=coord(:,1)
        if (kk.eq.2) kkk=2
        if (kk.eq.3) kkk=1
     end if
     do i=1,natom
        if ((charge(i).gt.1).and.(coord(i,1).gt.0.d0)) then  ! TEST
           !plus
           do k=1,2
              if (k.eq.1) signo=1!signo=1 -> add to coord of nucleus
              if (k.eq.2) signo=-1!signo=-1 -> subtract from coord of nucleus
              do j=1,f
                 if ((kk.eq.1).and.(tempa(i).gt.R_cutoff)) cycle!in r: if zone boundary bigger than r_cutoff, then it's not a nucleus, so don't need a k-shell radius

                 if ((signo*r(i)*(tempa(i)-tempb(j)).lt.0.2d0)&
                      &.and.(signo*r(i)*(tempa(i)-tempb(j)).ge.0.19d0) ) then

                    stat=.FALSE.!if there's already a zone boundary in the region, don't put another one there.
                    exit
                 else if ((j==f).and.(tempa(i)-signo*0.2d0/r(i).lt.Pi*kkk)) then!otherwise, put a k-shell node at the requisite boundary
                    tempb(f+1)=tempa(i)-signo*0.5d0/r(i)
                    !---------------------------------------------------
                    if ((tempa(i).eq.0.d0).and.(kk.eq.2).and.(tempb(f+1).lt.0.d0).and.(.not.control_zero_phi)) then
                       !special case: if phi coord close to zero, then put k-shell boundary at 2*Pi-value instead of negative value
                       tempb(f+1)=2.d0*Pi-abs(tempb(f+1))
                       control_zero_phi=.TRUE.
                    end if
                    !---------------------------------------------------
                    if (tempb(f+1).le.0.d0) then
                       tempb(f+1)=0.d0
                       stat=.FALSE.!ensure there are no negative values in the grid
                       exit
                    else
                       stat=.TRUE.
                    end if
                 else
                 end if
              end do
              if (stat) then
                 f=f+1
                 stat=.FALSE.
              end if
           end do
        end if
     end do
     if (kk.gt.1) then
        write(13,*)'tempb',kk,tempb(:)/Pi
     else
        write(13,*)'tempb',kk,tempb(:)
     end if
     f_store(kk)=f+1
     coord_a(:,kk)=tempb(:)
     !-------------------------------------------------
  end do


  !after setting up zone boundaries, want to sort boundaries in each dimension so that they lie in ascending order

  ! Ordering of coordinate vectors
  do kk=1,3
     do i=1,f_store(kk)
     write(13,*)'r',coord_a(i,kk),i,kk
     end do

     call sort_shell(coord_a(:,kk),dummy,f_store(kk))
     f_store(kk)=f_store(kk)-1
     do i=1,f_store(kk)
     write(13,*)'r_a',coord_a(i,kk),i,kk
     end do
  end do

  ! Construct cross-product basis
  f=(f_store(1)-1)*(f_store(2)-1)*(f_store(3)-1)
  write(6,*)(f_store(i),i=1,3)
  allocate(cross_product(f,6))
  write(16,*)coord_a(f_store(1),1)
  write(16,*)f
  do i=1,f_store(1)-1
     do k=1,f_store(2)-1
        do j=1,f_store(3)-1

           f=(i-1)*(f_store(2)-1)*(f_store(3)-1)+(k-1)*(f_store(3)-1)+j
           cross_product(f,1)=coord_a(i,1)
           cross_product(f,2)=coord_a(i+1,1)
           cross_product(f,3)=coord_a(k,2)/Pi
           cross_product(f,4)=coord_a(k+1,2)/Pi
           cross_product(f,5)=coord_a(j,3)/Pi
           cross_product(f,6)=coord_a(j+1,3)/Pi
           !delta=sector spacing
           !in areas close to nuclei, sector spacing = zone spacing
           !in areas far from nuclei, zone spacing > delta_r_max, delta_theta_max, or delta_phi_max, so set sector spacing
           !equal to appropriate variable delta_*_max
           delta(1)=(cross_product(f,2)-cross_product(f,1))
           delta(2)=(cross_product(f,4)-cross_product(f,3))
           delta(3)=(cross_product(f,6)-cross_product(f,5))
           if (delta(2).gt.delta_phi_max) delta(2)=delta_phi_max
           if (delta(3).gt.delta_theta_max) delta(3)=delta_theta_max
           if (delta(1).gt.delta_r_max) delta(1)=delta_r_max
           ! External region:
           if (cross_product(f,2).gt.(R0-epsilon_custom)) delta(1)=delta_external_radius
           write(16,900)(cross_product(f,n),n=1,6),delta(1),delta(2),delta(3)
        end do
     end do
  end do
  write(6,*)f
  deallocate(x,y,z,charge,coord,coord_a,tempa,tempb)
  deallocate(cross_product,delta)
900 format(9f13.8)
end program CM


subroutine sort_shell(eigen_vect,open_vect,n)
  implicit none
  integer n
  real*8,dimension(n),INTENT(INOUT) :: eigen_vect
  real*8,dimension(:),INTENT(INOUT) :: open_vect
  integer :: i,j,inc,l,max
  real*8 :: V
  !n=size(eigen_vect)
  inc=1
  max=n
  do
     inc=3*inc+1
     if (inc > n) exit
  end do
  do
     inc=inc/3
     do i=inc+1,n
        v=eigen_vect(i)
        j=i
        do
           if (eigen_vect(j-inc) <= v) exit
           eigen_vect(j)=eigen_vect(j-inc)
           j=j-inc
           if(j <= inc) exit
        end do
        eigen_vect(j)=v
     end do
     if (inc <= 1) exit
  end do
  do j=1, max
  end do
end subroutine sort_shell

