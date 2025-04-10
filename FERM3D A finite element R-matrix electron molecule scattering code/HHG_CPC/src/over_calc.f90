!!!*************************************************************
! 文件/File: over_calc.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: over_calc.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


module over
  use nrtype,only : dbl,i4b
  real(kind=dbl),dimension (:,:),allocatable :: open_mat_v,overlap_v,sovrapp_v
  integer(kind=i4b),dimension(:,:),allocatable :: open_mat_c,overlap_c,sovrapp_c
end module over

subroutine over_calc
  use nrtype, only : dbl,i4b
  use over
  use Brep
  use calc_func
  use control
  use open_information 
  use Matrices, only : open_mat=>gamma,overlap=>area_overlap
  implicit none
  real(kind=dbl),parameter :: a=0.0d0,b=1.0d0
  real(kind=dbl) :: ss,x_coord,y_coord,z_coord,mult1,mult2
  integer(kind=i4b),dimension(:),allocatable :: Patch_counter
  integer(kind=i4b) :: f,s,counter2,columnB,columnA,node,long_index
  write(6,*)'restore=',restore
  !------------------------------
  ! This subroutine constructs the numerical spherical harmonics by diagonalising the L^2
  ! operator on the boundary of the spherical R-matrix box, makes use of sparse
  ! solver SuperLU for the eigenvalue problem
  !------------------------------
  !------------------------------
  ! Variables
  ! overlap_v,open_mat_v = matrices with values of the overlap and del^2 for L^2 psi=l(l+1)OV psi
  ! bionda2_boundary = identifies wich functions are open on the boundary
  !            (global node numeration)  
  ! bounda_boundary = identifies wich functions are open on the boundary
  !            (local node numeration)

  !------------------------------
  !Allocation
  allocate(overlap_v(p_Patch,36))
  allocate(overlap_c(p_Patch,1))
  allocate(open_mat_v(p_Patch,36))
  allocate(open_mat_c(p_Patch,36))


  allocate(Patch_counter(p_Patch))
  allocate(bounda_boundary(n,64))
  allocate(bionda2_boundary(0:p_Patch))
  !  allocate(open_info(p_Patch,9))
  allocate(long_open(p_Patch,2))
  !----------------
  !Initialization
  index=0
  long_index=0
  do i=1,p_Patch
     Patch_counter(i)=0
     overlap_c(i,1) = 0 
     long_open(i,1)=0
     long_open(i,2)=0 
     do j=1,36
        open_mat_c(i,j) = 0
        open_mat_v(i,j) = 0.0d0 
        overlap_v(i,j) = 0.0d0
     end do
     bionda2_boundary(i) = 'none'
  end do
  do element=1,n
     do i=1,64
        bounda_boundary(element,i)='none'
     end do
  end do

  !-----------------------------------------
  ! DESCRIPTION
  ! Store in array long_open the local and element indices of the functions
  ! that are open at the surface 
  ! In practice just the ones that have the radial part equal to the 4th
  ! polynomial, the only one that is nonzero on the surface
  !-----------------------------------------

  do element=1,n
     do i=1,64
        if (boundary(nodes(element,int((i+7)/8))).eq.'up') then
           !$$$$$$$$$
           !if (theta_glo(nodes(element,int((i+7)/8))).le.(1.0d0-epsilon)) then
           !$$$$$$$$$

           if ((i.eq.33).or.(i.eq.35).or.(i.eq.36).or.(i.eq.39).or.(i.eq.41).or.(i.eq.43).or.(i.eq.44).or.&
                &(i.eq.47).or.(i.eq.49).or.(i.eq.51).or.(i.eq.52).or.(i.eq.55).or.(i.eq.57).or.(i.eq.59).or.(i.eq.60).or.(i.eq.63)) then 

              Patch_counter(PatchVertex(element,i))=Patch_counter(PatchVertex(element,i))+1
              bounda_boundary(element,i)='up_fun'
              !write(210,*)Patch_counter(PatchVertex(element,i))
              bionda2_boundary(PatchVertex(element,i))='up_fun'

              long_index=long_index+1
              long_open(long_index,1)=element
              long_open(long_index,2)=i

              if (Patch_counter(PatchVertex(element,i)).eq.1) then
                 index=index+1
                 !                 open_info(index,1)=PatchVertex(element,i)
                 !                 open_info(index,2)=element
                 !                 open_info(index,3)=i
                 !                 write(14,*)element,i,PatchVertex(element,i)
              else 
                 do j=1,index
                    !                    if (open_info(j,1).eq.PatchVertex(element,i)) then
                    select case (Patch_counter(PatchVertex(element,i)))
                    case(2)
                       !                          open_info(j,4)=element
                       !                          open_info(j,5)=i
                       !$                          write(15,*)element,i,PatchVertex(element,i)
                    case(3)
                       !                          open_info(j,6)=element
                       !                          open_info(j,7)=i
                       !$                          write(16,*)element,i,PatchVertex(element,i)
                    case(4)
                       !                          open_info(j,8)=element
                       !                          open_info(j,9)=i
                       !$                          write(17,*)element,i,PatchVertex(element,i)
                    end select
                    !                    end if
                 end do
              end if
           end if
           !$$$$$$$$
           !end if
           !$$$$$$$$

        end if
     end do
  end do
  max_index=index
  long_max=long_index
  !$write(119,*)long_index
  !do i=1,max_index
  !write(211,*)open_info(i,1),open_info(i,2),open_info(i,3),open_info(i,4),open_info(i,5),open_info(i,6),open_info(i,7),open_info(i,8),open_info(i,9)

  !end do

  write(6,*)'max_index',max_index,long_index
  !--------------------------------
  ! Construction of the overlap and del^2 matrices
  !  
  !--------------------------------
  do f=1,2
     if (f.eq.1) then
        choice_open='open'
     else if (f.eq.2) then 
        choice_open='overlap'
     end if
     if (f.eq.2) then
        do i=1,p_Patch
           do j=1,36
              open_mat_c(i,j)=0.0d0
           end do
        end do
        ! Loop over the elements
        do element=1,n

           call c_func

           ! Loops over the 64x64 functions inside the element
           do i=1,64
              do j=1,64
                 if ((bounda_boundary(element,i).eq.'up_fun').and.(bounda_boundary(element,j).eq.'up_fun')) then
                    call param_integr(i,k_ri,k_phii,k_thetai,mult1)
                    call param_integr(j,k_rj,k_phij,k_thetaj,mult2)
                    call overlap_2d(a,b,ss)
ss=ss*mult1*mult2
                    ! using a sparse format to store the matrices (the entries are always less or
                    ! equal to 36 per row)
                    do m=1,36
                       if (open_mat_c(PatchVertex(element,i),m).ne.0) then
                          if (open_mat_c(PatchVertex(element,i),m).eq.PatchVertex(element,j)) then 
                             overlap_v(PatchVertex(element,i),m)=overlap_v(PatchVertex(element,i),m)+ss
                             !                             write(340,*)overlap_c(PatchVertex(element,i),m)
                             exit
                          end if
                       else 
                          overlap_v(PatchVertex(element,i),m)=overlap_v(PatchVertex(element,i),m)+ss 
                          open_mat_c(PatchVertex(element,i),m)=PatchVertex(element,j)
                          exit
                       end if
                    end do
                 end if
              end do
           end do
        end do
        write(6,*)'done overlap'
     else
        do element=1,n
           call c_func
           do i=1,64
              do j=1,64
                 if((bounda_boundary(element,i).eq.'up_fun').and.(bounda_boundary(element,j).eq.'up_fun')) then
                    call param_integr(i,k_ri,k_phii,k_thetai,mult1)
                    call param_integr(j,k_rj,k_phij,k_thetaj,mult2)
                    call overlap_2d(a,b,ss)
ss=ss*mult1*mult2
                    do m=1,36
                       if (open_mat_c(PatchVertex(element,i),m).ne.0) then
                          if (open_mat_c(PatchVertex(element,i),m).eq.PatchVertex(element,j)) then
                             open_mat_v(PatchVertex(element,i),m)=open_mat_v(PatchVertex(element,i),m)+ss
                             !                             write(342,*)open_mat_c(PatchVertex(element,i),m),PatchVertex(element,i),m,element,i,j
                             exit
                          end if
                       else
                          open_mat_v(PatchVertex(element,i),m)=open_mat_v(PatchVertex(element,i),m)+ss
                          open_mat_c(PatchVertex(element,i),m)=PatchVertex(element,j)
                          !                          write(343,*)open_mat_c(PatchVertex(element,i),m),PatchVertex(element,i),m,element,i,j
                          exit
                       end if
                    end do
                 end if
              end do
           end do
        end do
     end if
  end do
  !----------------
  !*COMM* Checking 
  !                
  !  do i=1,p_Patch
  !     do j=1,36
  !        if (open_mat_c(i,j).ne.0) then
  !           write(216,900)overlap_v(i,j),i,open_mat_c(i,j)
  !        end if
  !     end do
  !  end do
  !900 format(f12.5,i6)
  !---------------
  !write(214,*)overlap_v(159,159)
  !------------------------
  !---------------------------------------------------------------------
  !This subroutine solves the motion on the circular boundary
  ! (and finds the surface harmonics)
  !---------------------------------------------------------------------

  !---------------
  !*COMM* Coordinates to plot the harmonics
  !                                        
  ! do i=1,p
  !     if (bionda2_boundary(8*i-7).eq.'up_fun') then
  !        x_coord=r_glo(i)*sin(theta_glo(i))*cos(phi_glo(i))
  !        y_coord=r_glo(i)*sin(theta_glo(i))*sin(phi_glo(i))
  !        z_coord=r_glo(i)*cos(theta_glo(i))
  !        write(3889,*)x_coord,y_coord,z_coord
  !     else 
  ! x_coord=r_glo(i)*sin(theta_glo(i))*cos(phi_glo(i))
  !        y_coord=r_glo(i)*sin(theta_glo(i))*sin(phi_glo(i))
  !        z_coord=r_glo(i)*cos(theta_glo(i))
  !        write(3890,*)x_coord,y_coord,z_coord
  !    end if
  !  end do  
  !---------------
  deallocate(Patch_counter)
  call sort(p_Patch)

  !-------------------------
  ! Call  sparse solvers to diagonalize
  !-------------------------

  call diag_sparse(p_Patch)
  deallocate(boundary)
end subroutine over_calc



! 2D Gauss-Legendre integration adapted


subroutine overlap_2d(x1,x2,ss)
  use nrtype, only : dbl
  implicit none
  real(kind=dbl) ss,x1,x2,overlap_h
  external overlap_h
  !This subroutine and the followings perform the 2D integrals (Numerical Recipes)

  call overlap_x(overlap_h,x1,x2,ss)
  return
end subroutine overlap_2d


function overlap_f(yy)
  use nrtype, only : dbl
  use control
  implicit none
  real(kind=dbl) overlap_f,yy,overlap_func,x,y
  common /xy/ x,y
  y=yy
  overlap_f=overlap_func(x,y)
  return
end function overlap_f



function overlap_h(xx)
  use nrtype, only : dbl
  implicit none
  real(kind=dbl) overlap_h,xx,overlap_f,y1,y2,x,y,ya,yb
  external overlap_f
  common /xy/ x,y
  real(kind=dbl) ss
  x=xx
  ya=0.0d0
  yb=1.0d0
  call overlap_y(overlap_f,ya,yb,ss)
  overlap_h=ss
  return
end function overlap_h





subroutine overlap_x(overlap_func,a,b,ss)
  use nrtype, only : dbl,i4b
  implicit none
  real(kind=dbl) a,b,ss,overlap_func
  external overlap_func
  integer(kind=i4b) j
  real(kind=dbl) dx,xm,xr,w(2),x(2)
  save w,x
  !  data w/.1279381955463162,.12583745654311546,.1216704731175888,.115505668233865,.10744427028407229,.09761865224963151,.0861901617116546,.07334648165006656,.05929858503369955,.044277439867617,.0285313860439130,.0123412297184586/

  !  data x/.06405689286260563,.1911188674736163,.3150426796961634,.4337935076260452,.5454214713888396,.6480936519369755,.7401241915785544,.820001985973903,.8864155270044011,.9382745520027328,.9747285559713095,.9951872199970213/

  ! data w/.295524224714753,.269266719309996,.219086362515982,.149451349150581,.066671344308688/
  ! data x/.148874338981631,.433339534129247,.679409568299024,.865063366688985,.973906528517172/
  data w/.652145154862546,.347854845137454/
  data x/.339981043584856,.861136311594053/
  xm=0.5d0*(b+a)
  xr=0.5d0*(b-a)
  ss=0.0d0
  do j=1,2
     dx=xr*x(j)
     ss=ss+w(j)*(overlap_func(xm+dx)+overlap_func(xm-dx))
  enddo
  ss=xr*ss
  return
end subroutine overlap_x



subroutine overlap_y(overlap_func,a,b,ss)
  use nrtype, only : dbl,i4b
  implicit none
  real(kind=dbl) a,b,ss,overlap_func
  external overlap_func
  integer(kind=i4b) j
  real(kind=dbl) dx,xm,xr,w(2),x(2)
  save w,x
  !  data w/.1279381955463162,.12583745654311546,.1216704731175888,.115505668233865,.10744427028407229,.09761865224963151,.0861901617116546,.07334648165006656,.05929858503369955,.044277439867617,.0285313860439130,.0123412297184586/

  !  data x/.06405689286260563,.1911188674736163,.3150426796961634,.4337935076260452,.5454214713888396,.6480936519369755,.7401241915785544,.820001985973903,.8864155270044011,.9382745520027328,.9747285559713095,.9951872199970213/
  !  data w/.295524224714753,.269266719309996,.219086362515982,.149451349150581,.066671344308688/
  !  data x/.148874338981631,.433395394129247,.679409568299024,.865063366688985,.973906528517172/
  data w/.652145154862546,.347854845137454/
  data x/.339981043584855,.861136311594053/
  xm=0.5d0*(b+a)
  xr=0.5d0*(b-a)
  ss=0.0d0
  do j=1,2
     dx=xr*x(j)
     ss=ss+w(j)*(overlap_func(xm+dx)+overlap_func(xm-dx))
  enddo
  ss=xr*ss
  return
end subroutine overlap_y





real(kind=8) function overlap_func(x,y)
  use nrtype, only : dbl
  use Calc_func
  use over
  use control,only : choice_open
  implicit none
  real(kind=8),intent(in) :: x,y
  if (choice_open.eq.'open') then
     !--------------
     !write(36,*)sin(angle_theta(x)),an_theta
     overlap_func=R0*R0*an_phi*an_theta*(sin(angle_theta(x))*der_theta(x,k_thetai)*der_theta(x,k_thetaj)*bas_phi(y,k_phii)*&
          &bas_phi(y,k_phij)/(an_theta*an_theta)+der_phi(y,k_phii)*der_phi(y,k_phij)*bas_theta(x,k_thetai)*bas_theta(x,k_thetaj)/(sin(angle_theta(x))*an_phi*an_phi))

     !     overlap_func=R0*R0*an_phi*an_theta*(sin(angle_theta(x))*der_thetai(x,i)*der_thetaj(x,j)*bas_phii(y,i)*bas_phij(y,j)/(an_theta*an_theta)+der_phii(y,i)*der_phij(y,j)*bas_thetai(x,i)*bas_thetaj(x,j)/(sin(angle_theta(x))*an_phi*an_phi))

  else
     overlap_func = (an_theta*an_phi)*(bas_theta(x,k_thetai)*bas_theta(x,k_thetaj)*bas_phi(y,k_phii)*bas_phi(y,k_phij)) *R0 *R0* sin(angle_theta(x))

     !     overlap_func = (an_theta*an_phi) *(bas_thetai(x,i)*bas_thetaj(x,j)*bas_phii(y,i)*bas_phij(y,j)) *R0 *R0* sin(angle_theta(x))

  end if
end function overlap_func











subroutine sort(p_Patch)
  use over
  use nrtype
  implicit none
  real(kind=dbl) :: prov,provola
  integer(kind=i4b) :: i,j,inc,n,max,s
  integer(kind=i4b),intent(in) :: p_Patch
  real(kind=dbl) :: V
  !-------------------------------
  ! Sorting subroutine to put matrix element in 2D problem in order for the
  ! reduced set of 2D nodes
  !-------------------------------
  n=36
  inc=1
  max=n
  do s=1,p_Patch
     do
        inc=3*inc+1
        if (inc > n) exit
     end do
     do
        inc=inc/3
        do i=inc+1,n
           if (open_mat_c(s,i).ne.0) then
              v=open_mat_c(s,i)
              prov=open_mat_v(s,i)
              provola=overlap_v(s,i)
              j=i
              do
                 if (open_mat_c(s,j-inc) <= v) exit
                 open_mat_c(s,j)=open_mat_c(s,j-inc)
                 open_mat_v(s,j)=open_mat_v(s,j-inc)
                 overlap_v(s,j)=overlap_v(s,j-inc)
                 j=j-inc
                 if(j <= inc) exit
              end do
              open_mat_c(s,j)=v
              open_mat_v(s,j)=prov
              overlap_v(s,j)=provola
           end if
        end do
        if (inc <= 1) exit
     end do
  end do
  do j=1, max
  end do
end subroutine sort
