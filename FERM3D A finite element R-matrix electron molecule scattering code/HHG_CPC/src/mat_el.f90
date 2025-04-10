!!!*************************************************************
! 文件/File: mat_el.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: mat_el.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


subroutine mat_el
  use nrtype,only : dp,Pi
  use Brep
  use Solve
  use Calc_func
  use Control, only : choice_partition,partition,choice_matrix,calculation_type,DFT,option_gauge
  use Matrices
  use Open_information
  use over,only : overlap_c
  use bsp_info
  use V_setup
  use gauss3dint
  !  use function_integr
  use function_integr_tab
  implicit none
  INTEGER, PARAMETER :: I44 = SELECTED_INT_KIND(9)
  integer(kind=i44) tempo_max
  integer(kind=i44) :: l,f,q,ind,k,kk,k_ind,kk_ind,qp
  real(kind=dp) :: basis,func
  real(kind=dp),parameter :: a=0.0d0,b=1.0d0
  real(kind=dp),dimension(3) :: aa,bb
  real(kind=dp) :: ss,r,phi,ss_corr
  real(kind=dp),parameter :: tolerance= 1.0d-14
  integer(kind=i44) :: nw,nEmax,start,symi,symj,ni,nj,lmax_loc,emme,mx,patch
  integer(kind=i44),parameter :: mmax=100,dime=216
  real(kind=dp) :: Emax,Emin,Delta_E
  logical,dimension(:,:),allocatable :: opench
  real(kind=dp), dimension(:,:),allocatable :: Ekin,T
  real(kind=dp), dimension(:), allocatable :: Energy
  integer(kind=i44),dimension(:), allocatable :: max_open_array,l_max
  real(kind=dp), dimension(:,:), allocatable ::area_overlapE,area_overlapU,area_overlap_open_closedE,&
       &area_overlap_open_closedU,area_overlap_open_openE,area_overlap_open_openU
  real(kind=dp),dimension(:),allocatable :: result_gamma,resultE,resultU
  integer(kind=i44),dimension(:),allocatable :: beta_v,singularity(:)
  integer(kind=i44) :: vect_ind,max_ind_sing
  integer(kind=i44) :: na,ll,cube_num_pts_sing,sing_sectors
  real(kind=dp) ::  r_glo_in,phi_glo_in,theta_glo_in,mult1,mult2
  real(kind=dp),parameter :: errass = 1.0d-8, errrel = 1.0d-6, dist=1.1763d0,E_last_bound=.57d0
  real(kind=dp) :: errest,separate_integral
  integer(kind=i44) :: max_num,num_partitions,square_num_func
  integer(kind=i44), parameter :: maxfcn = 60000,num_pts_sing=16,num_pts=4
  character(10) :: sing_elim_choice
  character(len=5), dimension(:),allocatable :: check_nuclei
  logical :: control_cycle,if_stack
  external ::r_glo_in,phi_glo_in,theta_glo_in,if_stack,separate_integral
  !---------------
  !Ionization Potentials
  ! H2 = .564
  ! N2 = .573
  ! CO2 = .50675
  !---------------

  !-----------------------------------------------------------------------
  ! DESCRIPTION
  !
  ! This subroutine calculates the matrix elements of the gamma matrix
  ! Variables: 
  ! l = index
  ! basis = function, defines the four basic polynomials
  ! func_* = functions,integrands of the 3D integrals
  ! a,b = extremes of the interval of integration (from 0 to 1 always)
  ! ss = result of the 3D integration within a sector
  ! nw = unit from where Energy data is read
  ! nE = index of number of steps in Energy
  ! nEmax = max number of steps in Energy
  ! mmax = max value of m allowed by the program
  ! E = value of the Energy
  ! Emax,Emin = max and min value of E, read from nw
  ! Delta_E = step in Energy
  ! opench = variable that says if a channel is open or closed
  ! Ekin = Kinetic Energy for each channel
  ! T = Threshhold for each channel
  !-----------------------------------------------------------------------

  write(6,*)'p_Patch',p_Patch
  nw=12
  max_num=num_pts**3
  open(unit=nw,file='input_control.dat',status='old',action='read')
  do i=1,9
  read(nw,*)
  end do
  read(nw,*)Emax
  read(nw,*)Emin
  read(nw,*)nEmax
  close(nw)

  ! Input checking
  if (Emax.lt.0.d0) then 
    write(6,*)'error, the maximum energy is negative'
    stop
  end if
  if (Emin.lt.0.d0) then
    write(6,*)'error, the minimum energy is negative'
    stop
  end if
  if (nEmax.lt.1) then
    write(6,*)'error, the number of  energy points is lower than 1'
    stop
  end if
  !----------------------------------------------------------------------------

  Delta_E=(Emax-Emin)/dble(nEmax)
  !read(nw,*)Delta_E 
  nEmax=anint((Emax-Emin)/Delta_E)
  !----------------------------------------------------------------------------
  !Allocation
  allocate(opench(mmax+1,nEmax))
  allocate(Ekin(mmax+1,nEmax))
  allocate(T(mmax+1,nEmax))
  allocate(Energy(nEmax))
  allocate(max_open_array(nEmax))
  allocate(l_max(nEmax))
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! WARNING
  ! The open channels should not exceed the number of the open functions
  !----------------------------------------------------------------------------
  !open(1300,form='unformatted')
  !open(1301,form='unformatted')
  !open(1302,form='unformatted')
  !-------------

  !Initialization
  do i=1,nEmax
     max_open_array(i)=0
     Energy(i)=0.0d0
  end do
  !-------------

  do nE=1,nEmax
     E=Emin+nE*Delta_E
     Energy(nE)=E 
     max_open=0 
     do li=0,mmax
        T(li+1,nE)=li*(li+1)/R0**2.0d0
        Ekin(li+1,nE)=E-T(li+1,nE)
        if (Ekin(li+1,nE) .gt. 0.0d0) then
           opench(li+1,nE)=.TRUE.
           max_open=max_open+1 
        else
           opench(li+1,nE)=.FALSE.
        end if
     end do
!!$     l_max(nE)=14
!max_open+0
!!$     max_open_array(nE)=196
!l_max(nE)*l_max(nE)
     l_max(nE)=max_open+1
     max_open_array(nE)=l_max(nE)*l_max(nE)
!---------------------
l_max(nE)=max_l
max_open_array(nE)=max_l**2
!---------------------
  end do
  deallocate(opench)
  deallocate(Ekin)
  deallocate(T)

  !--------------------------------------------------------------------------------
  ! VARIABLES
  ! Allocation of matrices for sparse algorhitm
  ! gamma,open_closed,open_open = matrices to store GAMMA matrix elements (see
  ! Greene et al. Rev. Mod. Phys. 68,1015 (pag.1043)) for closed-closed,
  ! open-closed,open-open partitions
  ! area_overlap,area_overlap_open_closed,area_overlap_open_open= same thing
  ! but for the overlap matrices
  ! gamma_c = matrix to store the numbering of the matrix elements from 1 to
  ! 216 for each row for the sparse format
  ! area_overlap_c = same thing for the overlap matrix
  ! interp_V, interp_density = interpolated values of the potential and density
  ! for every sector and integration point (except coulomb singularities)
  ! Interp_V_sing = 
  ! Interp_density_sing = interpolated values of the potential and density
  ! for every sector and integration point on the sectors that include
  ! singularities
  !--------------------------------------------------------------------------------

  allocate(gamma(p_Patch,dime))
  allocate(area_overlap_c(p_Patch))
  allocate(gamma_c(p_Patch,0:dime))
  allocate(area_overlapE(p_Patch,dime))
  allocate(area_overlapU(p_Patch,dime))
  allocate(interp_V(n,max_num))
  allocate(interp_density(n,max_num))
  allocate(singularity(n))
  !allocate(wf_ind(p_Patch))
  !-----------------------------------------------------------------------
  !Initialization 
  !-----------------------------------------------------------------------

  !$$---------------------------------------------------------------------


  area_overlap_c(:)=0
  gamma(:,:)=0.d0
  gamma_c(:,:)=0
  interp_v(:,:)=0.0d0
  interp_density(:,:)=0.0d0
  singularity(:)=0

  !$$---------------------------------------------------------------------
  do nE=1,nEmax
     max_open=max_open_array(nE)
     max_open_total=max_open_array(nEmax)
     write(6,*)'max open channels total=',max_open_array(nEmax),'max num harmonics in grid',max_index
     if (max_open_array(nEmax).gt.max_index) then
        write(6,*)'Error: number of open channels greater than number of points on surface'
        stop
     end if
     E=Energy(nE)
     lmax=l_max(nE)
     !-------------------------
     !lmax=3
     !max_open=9
     !-------------------------
     !$$---------------------------------------------------------------------
     write(6,*)'max_open channels ',max_open
     write(6,*)'at energy point nE=',nE
     allocate(area_overlap_open_closed(p_Patch,max_open))
     allocate(area_overlap_open_open(max_open,max_open))
     area_overlap_open_closed(:,:) = 0.0d0
     area_overlap_open_open(:,:)=0.0d0


     counter=0
     if (nE.eq.1) then
        tempo_max=max_open
        max_open=max_open_array(nEmax)
        allocate(open_closed(p_Patch,max_open))
        allocate(open_open(max_open,max_open))
        allocate(area_overlap_open_closedE(p_Patch,max_open))
        allocate(area_overlap_open_closedU(p_Patch,max_open))
        allocate(area_overlap_open_openE(max_open,max_open))
        allocate(area_overlap_open_openU(max_open,max_open))
        open_closed(:,:) = 0.0d0
        area_overlap_open_closedE(:,:) = 0.0d0
        area_overlap_open_closedU(:,:) = 0.0d0
        open_open(:,:)=0.0d0
        area_overlap_open_openE(:,:)=0.0d0
        area_overlap_open_openU(:,:)=0.0d0
        area_overlapE(:,:)=0.0d0
        area_overlapU(:,:)=0.0d0
        !call potential
        !XXXXXXXXXXXXXX
        if (calculation_type.eq.'bound') then
           num_partitions=1
        else 
           num_partitions=3
        end if

        !YYY
     end if
     !YYY

     !-----------------------
     ! **COMMENT
     ! Loop over partitions (closed-closed,open-closed,open-open)
     !-----------------------

     do partition=1,num_partitions
        !XXXXXXXXXXXXXX
        write(6,*)'partition',partition
        select case (partition)
        case(1)
           choice_partition='closed-closed'
           write(6,*)'closed-closed'  
        case(2) 
           choice_partition='open-closed'
           write(6,*)'open-closed' 
        case(3)
           choice_partition='open-open'
           write(6,*)'open-open'
        end select
        !YYYY 
        if (nE.eq.1) then
           !YYYY
           choice_func='func_en' 
           if (partition.eq.1) then
              choice_matrix='potential_empty'



              !*******************************************************************************************
              !--------------------------------------------------
              ! **COMM**
              ! Construction of indexing for sparse matrices
              !-------------------------------
              ! DESCRIPTION
              ! Usage of the sparse format to store matrix elements,
              ! row indices related to global indices with matrix gamma_c,
              ! gamma_c stores in advance the whole pattern of indices for the rows of
              ! gamma_v (closed-closed part)
              ! array area_overlap_c counts global index of function for only closed-closed functions
              !-------------------------------

              ind =0
              do j=1,p_Patch
                 if (bionda2_boundary(j).ne.'up_fun') then
                    ind=ind+1
                    area_overlap_c(j)=ind
                 end if
              end do
             ! Wavefunction indexing
             ! wf_ind(:)=0
             ! do i=1,ind
             !    do j=1,p_Patch
             !      if (area_overlap_c(j)==i) then
             !         wf_ind(i)=j
             !         write(35,*)wf_ind(i),i
             !      end if
             !    end do
             ! end do
             !-----------------------
              do element=1,n
                 do i=1,64
                    if (bounda_boundary(element,i).ne.'up_fun') then
                       do j=1,64
                          if (bounda_boundary(element,j).ne.'up_fun') then
                             f=0     
                             do q=1,dime
                                if (gamma_c(PatchVertex(element,i),q).ne.0) then
                                   if(gamma_c(PatchVertex(element,i),q).eq.PatchVertex(element,j)) then
                                      exit
                                   end if
                                else
                                   gamma_c(PatchVertex(element,i),q)=PatchVertex(element,j)
                                   f=f+1
                                   gamma_c(PatchVertex(element,i),0)=f 
                                   exit
                                end if
                             end do
                          end if
                       end do
                    end if
                 end do
              end do
           end if


           !---------------------------------------
           ! Sorting of the elements of gamma_c according to
           ! their global indices in the rows of the matrix
           !---------------------------------------

           call sort3(p_Patch)
           write(6,*)'restore=',restore

           !ZZZZ
           !ZZZZ
           !--------------------------------------------------
           !*********************************************************************************************
           if (choice_partition.eq.'closed-closed') then

              !-----------------------------
              ! Setting up different treatment for coulomb singularities (more integration
              ! points, different arrays where we store the needed V, density values)
              !-----------------------------

              !-----------------------------
              ! bas_tab = array to store values of the basis functions polynomials (to avoid
              ! recalculating them all the time) 
              ! der_tab = same thing with the derivatives
              ! interp_V_sing, interp_dens_sing = array to store values of V, density for
              ! each integration point in the sectors with singularities
              !-----------------------------

              cube_num_pts_sing=num_pts_sing**3
              sing_sectors=4*nAtom*phi_sectors
              write(6,*)"num Legendre points on singularity",cube_num_pts_sing,sing_sectors
              allocate(vect_sing(sing_sectors))
              allocate(interp_V_sing(sing_sectors,cube_num_pts_sing))
              allocate(interp_density_sing(sing_sectors,cube_num_pts_sing))
              allocate(val_array(max_num))
              allocate(val_array_sing(cube_num_pts_sing)) 
              allocate(check_nuclei(nAtom))
              check_nuclei(:)='no'
              do i=1,sing_sectors
                 vect_sing(i) = 0
                 do j=1,cube_num_pts_sing
                    interp_V_sing(i,j) = 0.0d0
                    interp_density_sing(i,j) = 0.0d0
                 end do
              end do
              ind_sing=0
              do element=1,n
                 num=0
                 if (element.eq.1) then
                    allocate(bas_tab(4,num_pts))
                    allocate(bas_tab_sing(4,num_pts_sing))
                    allocate(der_tab(4,num_pts))
                    allocate(der_tab_sing(4,num_pts_sing))
                    countx=0
                    ss=qgss3d(func_tab_sing,a,b,a,b,a,b,num_pts_sing)
                    countx=0
                    ss=qgss3d(func_tab,a,b,a,b,a,b,num_pts)
                    !---------------------------------------------------------
                    square_num_func=64**2
                    allocate(mem_tab(square_num_func,max_num))
                    count_func=0
                    do i=1,64
                       call param_integr(i,k_ri,k_phii,k_thetai,mult1)
                       do j=1,64
                          call param_integr(j,k_rj,k_phij,k_thetaj,mult2)
                          countx=0
                          count_func=count_func+1
                          count_pts=0
                          ss=qgss3d(func_mem_tab,a,b,a,b,a,b,num_pts)*mult1*mult2
                       end do
                    end do

                    !---------------------------------------------------------
                    if (DFT) then
                       call  grad_nabla_integr_prep(n,max_num,sing_sectors,cube_num_pts_sing)
                    end if


                 end if
                 num=0
                 if (mod(element,100).eq.0) write(6,*)'potential empty',element
                 call c_func
                 control_cycle=.FALSE.
                 !-------------------
                 ! Store quantities to optimize integrals
                 !-------------------

                 do l=1,8
                    do i=1,nAtom
                       if (if_stack(r_glo(nodes(element,l)),theta_glo(nodes(element,l)),phi_glo(nodes(element,l)),&
                            &epsilon_custom,position(i),angle(i),angle_position_phi(i))) then
                          singularity(element)=1
                          write(6,*)'inside',position(i),angle(i),angle_position_phi(i)
                          check_nuclei(i)='yes'
                          ind_sing=ind_sing+1
                          write(6,*)ind_sing
                          vect_sing(ind_sing)=element
                          control_integr='sing'
                          countx=0
                          ss=qgss3d(func_potential_empty_sing_tab,a,b,a,b,a,b,num_pts_sing)
                          control_cycle=.TRUE.
                          exit
                       else if ((l.eq.8).and.(i.eq.nAtom)) then
                          control_integr='outside' 
                          countx=0
                          ss=qgss3d(func_potential_empty_tab,a,b,a,b,a,b,num_pts)
                          control_cycle=.TRUE.
                          exit
                       end if
                    end do
                    if (control_cycle) exit
                 end do
              end do
              write(6,*)'Atoms with grid points at nucleus positions'
              do i=1,nAtom
               if (position(i).eq.0.d0) check_nuclei(i)='yes'
               write(6,*)'atom',i,check_nuclei(i)
              end do
              !XXXXXXXXXX
              !*COMM* Deallocation of spline coefficients arrays for V and density
              if (option_gauge.ne.'acceleration') then
                deallocate(bscoef_V)
                deallocate(bscoef_density)
                if (DFT) then
                   deallocate(bscoef_grad)
                   deallocate(bscoef_nabla)
                end if
              end if
              !XXXXXXXXXX
           end if
           write(6,*)'after potential empty'
           max_ind_sing=ind_sing
           !***********************************************************************************************
           if (restore.eq.'no') then
              !---------------------------
              ! Construction of the gamma matrices
              !---------------------------
              !--------------------------------------------------
              choice_matrix='gamma'
              if (choice_partition.eq.'closed-closed') then
                 do element=1,n
                    !----------------------------------------
                    !$$$ Singularity elimination for 1/sin(theta) in phi term of kinetica energy 
                    !$$$ (does not work yet)
                    !----------------------------------------
                    !theta_in=Pi/2.d0
                    !sing_elim_choice='func'
                    !do l=1,8
                    !if ((theta_glo(nodes(element,l)).le.epsilon_custom).or.(theta_glo(nodes(element,l)).ge.(Pi-epsilon_custom))) then
                    !theta_in=theta_glo(nodes(element,l))
                    !sing_elim_choice='func_sing'
                    !exit
                    !end if
                    !end do
                    !-----------------------------------------
                    if (mod(element,100).eq.0) write(6,*)'closed-closed gamma'
                    call c_func
                    do i=1,64
                       if (bounda_boundary(element,i).ne.'up_fun') then 
                          call param_integr(i,k_ri,k_phii,k_thetai,mult1)
                          do j=1,64
                             if (bounda_boundary(element,j).ne.'up_fun') then
                                call param_integr(j,k_rj,k_phij,k_thetaj,mult2)
                                do q=1,dime
                                   if(gamma_c(PatchVertex(element,i),q).eq.PatchVertex(element,j)) then
                                      !if (PatchVertex(element,i).ge.gamma_c(PatchVertex(element,i),q)) then
                                      countx=0
                                      !-----------------------------
                                      !$$$  Phi singularity elimination
                                      !sing_elim_choice='func'
                                      !ss=separate_integral(a,b,8,element,theta_in,sing_elim_choice)
                                      !-----------------------------
                                      !-----------------------------
                                      !$$$  Separate form of the integrals
                                      ss=separate_integral(a,b,4,element)*mult1*mult2
                                      !$$$  3D integration (slow)
                                      !ss=qgss3d(func_closed_closed_gamma_tab,a,b,a,b,a,b,num_pts)
                                      !-----------------------------
                                      gamma(PatchVertex(element,i),q)=gamma(PatchVertex(element,i),q)+ss
                                      !end if
                                      exit
                                      !end if
                                   end if
                                end do
                             end if
                          end do
                       end if
                    end do
                 end do
                 !**********************************************************************************
              else if (choice_partition.eq.'open-closed') then 
                 do element=1,n  
                    if (mod(element,100).eq.0)write(6,*)'open-closed gamma'
                    call c_func 
                    do i=1,64  
                       if (bounda_boundary(element,i).ne.'up_fun') then 
                          call param_integr(i,k_ri,k_phii,k_thetai,mult1)
                          do index=1,long_max
                             if (element.eq.(long_open(index,1))) then
                                beta=long_open(index,2)
                                call param_integr(beta,k_rj,k_phij,k_thetaj,mult2)
                                countx=0
                                ss=separate_integral(a,b,4,element)*mult1*mult2
                                !ss=qgss3d(func_open_closed_gamma_tab,a,b,a,b,a,b,num_pts)
                                do m=1,max_open
                                   open_closed(PatchVertex(element,i),m)=open_closed(PatchVertex(element,i),m)+&
                                        &ss*open_vect(overlap_c(PatchVertex(element,beta),1),m)
                                end do
                             end if
                          end do
                       end if
                    end do
                 end do
                 !*********************************************************************************************
              else  if (choice_partition.eq.'open-open') then

                 !*COMM* Acceleration of the open_open part of the matrices
                 !*COMM* Store the harmonics in an array 

                 !-----------------------------------------------------
                 ! VARIABLES
                 ! alpha_v = number of the first function (local numeration in the sector)
                 ! beta_v = number of the second function (local numeration in the sector)
                 ! countx = counter of x points for integration
                 ! rad_elem_v = element in which the open type functions are located
                 !-----------------------------------------------------
                 vect_ind=0
                 do index=1,long_max
                    do indexx=1,long_max
                       if(long_open(index,1).eq.long_open(indexx,1)) then
                          vect_ind=vect_ind+1
                       end if
                    end do
                 end do
                 vect_max=vect_ind
                 vect_ind=0
                 allocate(alpha_v(vect_max))
                 allocate(beta_v(vect_max))
                 allocate(rad_elem_v(vect_max))
                 allocate(result_gamma(vect_max))
                 do index=1,long_max
                    do indexx=1,long_max
                       if(long_open(index,1).eq.long_open(indexx,1)) then
                          vect_ind=vect_ind+1
                          alpha_v(vect_ind)=long_open(index,2)
                          beta_v(vect_ind)=long_open(indexx,2)
                          rad_elem_v(vect_ind)=long_open(index,1)
                          rad_elem=rad_elem_v(vect_ind)
                          alpha=alpha_v(vect_ind)
                          beta=beta_v(vect_ind)
                          call c_func
                          call param_integr(alpha,k_ri,k_phii,k_thetai,mult1)
                          call param_integr(beta,k_rj,k_phij,k_thetaj,mult2)
                          countx=0
                          ss=separate_integral(a,b,4,rad_elem)*mult1*mult2
                          !ss=qgss3d(func_open_open_gamma_tab,a,b,a,b,a,b,num_pts)
                          result_gamma(vect_ind)=ss
                       end if
                    end do
                 end do
!call natural_harmonics(long_open,max_open,long_max,max_index,vect_max,alpha_v,rad_elem_v)
                 !----------------------------------------------------------------------------
                 ! VARIABLES
                 ! open_vect = array of the coefficients of the spherical harmonics on the nodes
                 ! on the surface
                 ! open-open mtrix element is calculated storing the array of integrals
                 ! result_gamma(vect_ind) on all the points on the surface, then getting the
                 ! correct value by multipliyng by the coefficients of the spherical harmonics on
                 ! the same points (below)
                 !----------------------------------------------------------------------------

                 do m=1,max_open
                    if (mod(m,100).eq.1)write(6,*)'open-open gamma'
                    do mm=1,max_open
                       do vect_ind=1,vect_max
                          open_open(m,mm)=open_open(m,mm)+result_gamma(vect_ind)*&
                               & open_vect(overlap_c(PatchVertex(rad_elem_v(vect_ind),beta_v(vect_ind)),1),mm)*&
                               & open_vect(overlap_c(PatchVertex(rad_elem_v(vect_ind),alpha_v(vect_ind)),1),m)
                       end do
                    end do
                 end do
                 deallocate(result_gamma)
              end if
              !-----------------------------------------------------------------------
              !  Calculate Area Overlap matrix (E*psi*psi)
              !  Since the variation of the hara exchange potential with E has
              !  been expanded in Taylor series up to 1st term, the linear part
              !  (d(v_hara)/dE)
              !  is stored in this matrix too
              !  V_hara=V_hara(E=0)+E*V_hara'  
              !-----------------------------------------------------------------------

              ss=0    
              choice_matrix='overlap'
              counter=0
              !**************************************************************************************************** 
              if (choice_partition.eq.'closed-closed') then
                 Energy_Hara=2.0d0*(E+E_last_bound)
                 num=0
                 !--------------------------------
                 do element=1,n 
                    if (mod(element,100).eq.0)write(6,*)'closed-closed overlap',element
                    call c_func
                    do i=1,64
                       if (bounda_boundary(element,i).ne.'up_fun') then            
                          call param_integr(i,k_ri,k_phii,k_thetai,mult1)
                          do j=1,64
                             if (bounda_boundary(element,j).ne.'up_fun') then       
                                counter=counter+1
                                call param_integr(j,k_rj,k_phij,k_thetaj,mult2)
                                do q=1,dime 
                                   num=0
                                   if(gamma_c(PatchVertex(element,i),q).eq.PatchVertex(element,j)) then
                                      control_cycle=.FALSE.
                                      if (singularity(element)==1) then
                                         do ll=1,max_ind_sing
                                            if (vect_sing(ll).eq.element) then
                                               ind_sing=ll
                                            end if
                                         end do
                                         count_func=64*(i-1)+j
                                         count_pts=0
                                         countx=0
                                         ss=qgss3d(func_sing_0_tab,a,b,a,b,a,b,num_pts_sing)*mult1*mult2
                                         area_overlapE(PatchVertex(element,i),q)=area_overlapE(PatchVertex(element,i),q)-2*ss
                                         control_cycle=.TRUE.
                                         exit
                                      else 
                                         countx=0
                                         count_func=64*(i-1)+j
                                         count_pts=0
                                         ss= qgss3d(func_closed_closed_overlap_tab,a,b,a,b,a,b,num_pts)*mult1*mult2
                                         area_overlapE(PatchVertex(element,i),q)=area_overlapE(PatchVertex(element,i),q)-2*ss
                                         control_cycle=.TRUE.
                                      end if
                                      exit 
                                   end if
                                end do
                             end if
                          end do
                       end if
                    end do
                 end do
                 !********************************************************************************
              else  if (choice_partition.eq.'open-closed') then
                 do element=1,n  
                    if (mod(element,100).eq.0)write(6,*)'open-closed overlap'
                    call c_func
                    count_func=0
                    do i=1,64
                       if (bounda_boundary(element,i).ne.'up_fun') then
                          call param_integr(i,k_ri,k_phii,k_thetai,mult1)
                          do index=1,long_max
                             count_func=count_func+1
                             if (element.eq.(long_open(index,1))) then
                                num=0
                                beta=long_open(index,2)
                                call param_integr(beta,k_rj,k_phij,k_thetaj,mult2)
                                countx=0
                                ss=qgss3d(func_open_closed_overlap_tab,a,b,a,b,a,b,num_pts)*mult1*mult2
                                do m=1,max_open
                                   area_overlap_open_closedE(PatchVertex(element,i),m)=&
                                        & area_overlap_open_closedE(PatchVertex(element,i),m)&
                                        & -2*ss*open_vect(overlap_c(PatchVertex(element,beta),1),m)
                                end do
                             end if
                          end do
                       end if
                    end do
                 end do
                 !****************************************************************************
              else if (choice_partition.eq.'open-open') then
                 counter=0
                 ss=0.0d0

                 allocate(resultE(vect_max))
                 vect_ind=0
                 count_func=0
                 do index=1,long_max
                    do indexx=1,long_max
                       count_func=count_func+1
                       if(long_open(index,1).eq.long_open(indexx,1)) then
                          vect_ind=vect_ind+1
                          alpha_v(vect_ind)=long_open(index,2)
                          beta_v(vect_ind)=long_open(indexx,2)
                          rad_elem_v(vect_ind)=long_open(index,1)
                          rad_elem=rad_elem_v(vect_ind)
                          alpha=alpha_v(vect_ind)
                          beta=beta_v(vect_ind)
                          element=rad_elem
                          num=0
                          call c_func
                          call param_integr(alpha,k_ri,k_phii,k_thetai,mult1)
                          call param_integr(beta,k_rj,k_phij,k_thetaj,mult2)
                          countx=0
                          ss=qgss3d(func_open_closed_overlap_tab,a,b,a,b,a,b,num_pts)*mult1*mult2
                          resultE(vect_ind)=-2*ss
                       end if
                    end do
                 end do
                 do m=1,max_open
                    if (mod(m,100).eq.1)write(6,*)'open-open overlap'
                    do mm=1,max_open
                       do vect_ind=1,vect_max
                          area_overlap_open_openE(m,mm)=area_overlap_open_openE(m,mm)+resultE(vect_ind)*&
                               & open_vect(overlap_c(PatchVertex(rad_elem_v(vect_ind),beta_v(vect_ind)),1),mm)&
                               & *open_vect(overlap_c(PatchVertex(rad_elem_v(vect_ind),alpha_v(vect_ind)),1),m)
                       end do
                    end do
                 end do
                 deallocate(resultE)
              end if
              !-------------------
              !*COMM* Restore Part
              do i=1,p_Patch
                 do j=1,dime
                                        !write(1000,*)gamma(i,j),area_overlapE(i,j)
                 end do
                 do j=1,max_open
                                        !write(2000,*)open_closed(i,j),area_overlap_open_closedE(i,j)
                 end do
              end do
              do i=1,max_open
                 do j=1,max_open
                                        !write(3000,*)open_open(i,j),area_overlap_open_openE(i,j)
                 end do
              end do
              !--------------------
           end if
           !YYYY 
        end if
        !YYYY
        !--------------------------------------------
        !--------------------------
        ! Calculate the potential matrix
        !--------------------------
        if ((nE.eq.1).or.(exchange.ne.'slater')) then
           do i=1,3
              aa(i)=0.0d0
              bb(i)=1.0d0
           end do
           ss=0    
           choice_matrix='potential'
           counter=0
           if (choice_partition.eq.'closed-closed') then
              area_overlapU=0.0d0
              Energy_Hara=2.0d0*(E+E_last_bound)
              !rewind(1300)
              do element=1,n 
                 if (mod(element,100).eq.0)write(6,*)'closed-closed potential',element
                 call c_func
                 do i=1,64
                    if (bounda_boundary(element,i).ne.'up_fun') then            
                       call param_integr(i,k_ri,k_phii,k_thetai,mult1)
                       do j=1,64
                          if (bounda_boundary(element,j).ne.'up_fun') then       
                             counter=counter+1
                             call param_integr(j,k_rj,k_phij,k_thetaj,mult2)
                             do q=1,dime 
                                num=0
                                if(gamma_c(PatchVertex(element,i),q).eq.PatchVertex(element,j)) then
                                   control_cycle=.FALSE.
                                   if (singularity(element)==1) then
                                      do ll=1,max_ind_sing
                                         if (vect_sing(ll).eq.element) then
                                            ind_sing=ll
                                         end if
                                      end do
                                      count_func=64*(i-1)+j
                                      count_pts=0
                                      countx=0
                                      ss=qgss3d(func_sing_1_tab,a,b,a,b,a,b,num_pts_sing)*mult1*mult2
                                      area_overlapU(PatchVertex(element,i),q)=area_overlapU(PatchVertex(element,i),q)+2*ss
                                      control_cycle=.TRUE.
                                      exit
                                   else 
                                      countx=0
                                      count_func=64*(i-1)+j
                                      count_pts=0
                                      ss= qgss3d(func_closed_closed_potential_tab,a,b,a,b,a,b,num_pts)*mult1*mult2
                                      area_overlapU(PatchVertex(element,i),q)=area_overlapU(PatchVertex(element,i),q)+2*ss
                                      control_cycle=.TRUE.
                                   end if
                                   exit 
                                end if
                             end do
                          end if
                       end do
                    end if
                 end do
              end do
           else  if (choice_partition.eq.'open-closed') then
              area_overlap_open_closedU=0.0d0
              do element=1,n  
                 if (mod(element,100).eq.0)write(6,*)'open-closed potential'
                 !rewind(1301)
                 call c_func
                 count_func=0
                 do i=1,64
                    if (bounda_boundary(element,i).ne.'up_fun') then
                       call param_integr(i,k_ri,k_phii,k_thetai,mult1)
                       do index=1,long_max
                          count_func=count_func+1
                          if (element.eq.(long_open(index,1))) then  
                             num=0
                             beta=long_open(index,2)
                             call param_integr(beta,k_rj,k_phij,k_thetaj,mult2)
                             countx=0
                             ss=qgss3d(func_open_closed_potential_tab,a,b,a,b,a,b,num_pts)*mult1*mult2
                             do m=1,max_open
                                area_overlap_open_closedU(PatchVertex(element,i),m)=&
                                     & area_overlap_open_closedU(PatchVertex(element,i),m)+&
                                     & 2*ss*open_vect(overlap_c(PatchVertex(element,beta),1),m)
                             end do
                          end if
                       end do
                    end if
                 end do
              end do
              !*******************************************************************************
           else if (choice_partition.eq.'open-open') then
              area_overlap_open_openU=0.0d0 
              !rewind(1302)
              counter=0
              ss=0.0d0

              vect_ind=0
              do index=1,long_max
                 do indexx=1,long_max
                    if(long_open(index,1).eq.long_open(indexx,1)) then
                       vect_ind=vect_ind+1
                    end if
                 end do
              end do
              vect_max=vect_ind
              vect_ind=0
              allocate(resultU(vect_max))
              vect_ind=0
              count_func=0
              do index=1,long_max
                 do indexx=1,long_max
                    count_func=count_func+1
                    if(long_open(index,1).eq.long_open(indexx,1)) then
                       vect_ind=vect_ind+1
                       alpha_v(vect_ind)=long_open(index,2)
                       beta_v(vect_ind)=long_open(indexx,2)
                       rad_elem_v(vect_ind)=long_open(index,1)
                       num=0
                       rad_elem=rad_elem_v(vect_ind)
                       element=rad_elem
                       alpha=alpha_v(vect_ind)
                       beta=beta_v(vect_ind)
                       call c_func
                       call param_integr(alpha,k_ri,k_phii,k_thetai,mult1)
                       call param_integr(beta,k_rj,k_phij,k_thetaj,mult2)
                       countx=0
                       ss=qgss3d(func_open_closed_potential_tab,a,b,a,b,a,b,num_pts)*mult1*mult2
                       resultU(vect_ind)=2*ss
                    end if
                 end do
              end do
              do m=1,max_open
                 if (mod(m,100).eq.1)write(6,*)'open-open potential'
                 do mm=1,max_open
                    do vect_ind=1,vect_max
                       area_overlap_open_openU(m,mm)=area_overlap_open_openU(m,mm)+resultU(vect_ind)*&
                            & open_vect(overlap_c(PatchVertex(rad_elem_v(vect_ind),beta_v(vect_ind)),1),mm)*&
                            & open_vect(overlap_c(PatchVertex(rad_elem_v(vect_ind),alpha_v(vect_ind)),1),m)
                    end do
                 end do
              end do
              deallocate(resultU)
              !XXXXXXXXXXXXXXXXXXXXXX
           end if



           !--------------------------------------------
        end if
     end do
     if (nE.eq.1) then
        ! DEBUG
        !write(42,*)area_overlapE
        !write(43,*)area_overlap_open_closedE
        !write(44,*)area_overlap_open_openE
        !write(54,*)area_overlapU
        !write(55,*)area_overlap_open_closedU
        !write(56,*)area_overlap_open_openU
        !-----------------------------------
        max_open=tempo_max
        deallocate(singularity)
        !deallocate(alpha_v)
        deallocate(beta_v)
        !deallocate(rad_elem_v)
        ! New deallocations
        !deallocate(r_glo)
        !deallocate(theta_glo)
        !deallocate(nodes)
        deallocate(r_glo_in_tab)
        deallocate(theta_glo_in_tab)
        !deallocate(delel_r,delel_theta,delel_phi)
        call integral_release_memory(max_num,sing_sectors,cube_num_pts_sing)
        deallocate(bounda_boundary)
        deallocate(position,angle,angle_position_phi)
     end if

     !---------------------------------------------------------------
     ! Infinite wall potential on the boundary


     !---------------------------------------------------------------
     !------------------------------------
     ! COMMENT
     ! Summing matrices for solution of R-matrix  eigenvalue problem
     ! Summing matrices differently for bound state problem (only closed-closed
     ! partitions) or R-matrix scattering problem 
     !------------------------------------
     allocate(area_overlap(p_Patch,dime))
     area_overlap=0.0d0
     nnz=0
     do i=1,p_Patch
        do j=1,dime
           if (restore.eq.'yes') then
              read(1000,*)gamma(i,j),area_overlapE(i,j)
           end if
           if (calculation_type.eq.'bound') then
              gamma(i,j)=gamma(i,j)+area_overlapU(i,j)
              area_overlap(i,j)=-area_overlapE(i,j)
           else
              area_overlap(i,j)=area_overlapU(i,j)+E*area_overlapE(i,j) 
           end if
           if (((dabs(gamma(i,j)).ge.tolerance).or.(dabs(area_overlap(i,j)).ge.tolerance))) then
              nnz=nnz+1
           end if
        end do
     end do
     !call struct_mat
901  format(f12.5,2i6)
900  format(1f13.8)
     do i=1,p_Patch
        do m=1,max_open
           if (restore.eq.'yes') then
              read(2000,*)open_closed(i,j),area_overlap_open_closedE(i,j)
           end if
           area_overlap_open_closed(i,m)=area_overlap_open_closedU(i,m)+E*area_overlap_open_closedE(i,m)

           !write(9000+nE,900)open_closed(i,m), i,m
           !write(10000+nE,900)area_overlap_open_closed(i,m)
           !write(32,*)area_overlap_open_closed(i,m)
        end do
     end do
     do m=1,max_open
        do mm=1,max_open
           if (restore.eq.'yes') then
              read(3000,*)open_open(i,j),area_overlap_open_openE(i,j)
           end if
           area_overlap_open_open(m,mm)=area_overlap_open_openU(m,mm)+E*area_overlap_open_openE(m,mm)
           !write(5000+nE,900)open_open(m,mm)
           !write(6000+nE,900)area_overlap_open_open(m,mm) 
           !write(35,*)area_overlap_open_open(m,mm)
        end do
     end do


     select case (problem) 
        ! COMMENT Call the module procedure Solve to perform the linear algebra of the problem 
     case (1)
        !        call lin(gamma)

     case(2)
        ! Calls the subroutine that performs all the linear algebra
        ! for the eiganvalue problem solution

        call lin 

     end select
     !$$---------------------------------------------------------------------
     !call bspsync()
  end do

  ! Deallocations
  deallocate(l_max)
  deallocate(max_open_array)
  deallocate(Energy)
  deallocate(gamma,gamma_c,area_overlap_c,area_overlapE,area_overlapU)
  deallocate(area_overlap_open_closedE,area_overlap_open_openE)
  deallocate(area_overlap_open_closedU,area_overlap_open_openU)
  !deallocate(PatchVertex)
  deallocate(coord_x,coord_y,coord_z)
  !deallocate(bionda2_boundary,long_open)
  deallocate(xknot,yknot,zknot)
 !$$---------------------------------------------------------------------

end subroutine mat_el



!----------------------------------
! Functions to calculate finite elements polynomials (and derivatives)
!----------------------------------

real(kind=8) function  bas_r(x,m)
  use nrtype,only : dbl,i4b
  implicit none
  real(kind=dbl) basis
  real(kind=dbl) x2,x3
  real(kind=dbl),intent(in) ::  x
  integer(kind=i4b),intent(in) ::  m
  x2=x*x
  x3=x2*x
  if (m.eq.1) then
     bas_r=1.0d0-3.0d0*x2+2.0d0*x3
  else if (m.eq.2) then
     bas_r=x3-2.0d0*x2+x
  else if (m.eq.3) then
     bas_r=3.0d0*x2-2*x3
  else if (m.eq.4) then
     bas_r=x3-x2
  end if
end function bas_r


real(kind=8) function  bas_phi(x,m)
  use nrtype,only : dbl,i4b
  use calc_func,only : element
  use Brep,only : epsilon_custom
  implicit none
  real(kind=dbl) basis
  real(kind=dbl) x2,x3
  real(kind=dbl),intent(in) ::  x
  real(kind=dbl) theta_glo_in
  external theta_glo_in
  integer(kind=i4b),intent(in) ::  m
  x2=x*x
  x3=x2*x
  if (m.eq.1) then
     bas_phi=1.0d0-3.0d0*x2+2.0d0*x3
  else if (m.eq.2) then
     bas_phi=x3-2.0d0*x2+x
  else if (m.eq.3) then
     bas_phi=3.0d0*x2-2*x3
  else if (m.eq.4) then
     bas_phi=x3-x2
  end if
end function bas_phi

real(kind=8) function  bas_theta(x,m)
  use nrtype,only : dbl,i4b
  implicit none
  real(kind=dbl) basis
  real(kind=dbl) x2,x3
  real(kind=dbl),intent(in) ::  x
  integer(kind=i4b),intent(in) :: m
  x2=x*x
  x3=x2*x
  if (m.eq.1) then
     bas_theta=1.0d0-3.0d0*x2+2.0d0*x3
  else if (m.eq.2) then
     bas_theta=x3-2.0d0*x2+x
  else if (m.eq.3) then
     bas_theta=3.0d0*x2-2*x3
  else if (m.eq.4) then
     bas_theta=x3-x2
  end if
end function bas_theta






real(kind=8) function  der_r(x,m)
  use nrtype,only : dbl,i4b
  implicit none
  real(kind=dbl) basis 
  real(kind=dbl) der_bas,x2,x3
  real(kind=dbl),intent(in) ::  x
  integer(kind=i4b),intent(in) ::  m
  x2=x*x
  x3=x2*x
  if (m.eq.1) then
     der_r=-6.0d0*x+6.0d0*x2
  else if (m.eq.2) then
     der_r=3.0d0*x2-4.0d0*x+1.0d0
  else if (m.eq.3) then
     der_r=6.0d0*x-6.0d0*x2
  else if (m.eq.4) then
     der_r=3.0d0*x2-2.0d0*x
  end if
end function der_r


real(kind=8) function  der_phi(x,m)
  use nrtype,only : dbl,i4b
  use calc_func,only : element
  use Brep,only : epsilon_custom
  implicit none 
  real(kind=dbl) der_bas
  real(kind=dbl) x2,x3
  real(kind=dbl),intent(in) ::  x
  real(kind=dbl) theta_glo_in
  external theta_glo_in
  integer(kind=i4b),intent(in) ::  m
  x2=x*x
  x3=x2*x
  if (m.eq.1) then
     der_phi=-6.0d0*x+6.0d0*x2
  else if (m.eq.2) then
     der_phi=3.0d0*x2-4.0d0*x+1.0d0
  else if (m.eq.3) then
     der_phi=6.0d0*x-6.0d0*x2
  else if (m.eq.4) then
     der_phi=3.0d0*x2-2.0d0*x
  end if
end function der_phi


real(kind=8) function  der_theta(x,m)
  use nrtype,only : dbl,i4b
  implicit none 
  real(kind=dbl) der_bas
  real(kind=dbl) x2,x3
  real(kind=dbl),intent(in) :: x
  integer(kind=i4b),intent(in) ::  m
  x2=x*x
  x3=x2*x
  if (m.eq.1) then
     der_theta=-6.0d0*x+6.0d0*x2
  else if (m.eq.2) then
     der_theta=3.0d0*x2-4.0d0*x+1.0d0
  else if (m.eq.3) then
     der_theta=6.0d0*x-6.0d0*x2
  else if (m.eq.4) then
     der_theta=3.0d0*x2-2.0d0*x
  end if
end function der_theta





real(kind=8) function der_bas(x,i)
  use nrtype,only : dbl,i4b
  implicit none
  real(kind=dbl) x2
  real(kind=dbl),intent(in) ::  x
  integer(kind=i4b),intent(in) ::  i
  x2=x*x
  select case (i)
  case(1)
     der_bas=-6.0d0*x+6.0d0*x2
  case(2)
     der_bas=3.0d0*x2-4.0d0*x+1.0d0
  case(3)
     der_bas=6.0d0*x-6.0d0*x2
  case(4)
     der_bas=3.0d0*x2-2.0d0*x
  end select

end function der_bas


real(kind=8) function sec_derbas(x,i)
  use nrtype,only : dbl,i4b
  implicit none
  real(kind=dbl),intent(in) ::  x
  integer(kind=i4b),intent(in) ::  i
  select case (i)
  case(1)
     sec_derbas=-6.0d0+12.0d0*x
  case(2)
     sec_derbas=6.0d0*x-4.0d0
  case(3)
     sec_derbas=6.0d0-12.0d0*x
  case(4)
     sec_derbas=6.0d0*x-2.0d0
  end select

end function sec_derbas


!------------------
! Functions to calculate the values of radius, phi, theta for each integration
!  point
!------------------

real(kind=8) function radius(x)
  use nrtype,only : dbl,i4b
  use Brep
  use Calc_func, only : an_r,element,c_func
  use control
  use open_information
  implicit none
  real(kind=dbl),intent(in) :: x
  real(kind=dbl) :: r_glo_in
  external r_glo_in
  !call c_func
  if (choice_partition.ne.'open-open') then
     radius=r_glo_in(element)+x*an_r
  else
     radius=r_glo_in(rad_elem)+x*an_r
  end if
end function radius


real(kind=8) function angle_theta(y)
  use nrtype,only : dbl
  use Brep
  use Calc_func,only :an_theta,element,c_func
  use control
  use open_information
  implicit none
  real(kind=dbl),intent(in) :: y
  real(kind=dbl) :: theta_glo_in
  external theta_glo_in
  !call c_func
  if (choice_partition.ne.'open-open') then
     angle_theta=theta_glo_in(element)+y*an_theta
  else
     angle_theta=theta_glo_in(rad_elem)+y*an_theta
  end if
end function angle_theta

real(kind=8) function angle_phi(y)
  use nrtype,only : dbl
  use Brep
  use Calc_func,only :an_phi,element,c_func
  use control
  use open_information
  implicit none
  real(kind=dbl),intent(in) :: y
  real(kind=dbl) :: phi_glo_in
  external phi_glo_in
  call c_func
  if (choice_partition.ne.'open-open') then
     angle_phi=phi_glo_in(element)+y*an_phi
  else
     angle_phi=phi_glo_in(rad_elem)+y*an_phi
  end if
end function angle_phi



subroutine param_integr(m,k_r,k_phi,k_theta,mult)
  !*COMM* subroutine that chooses the number of the functions in the integral
  !*COMM* from the indices of the nodes 
use nrtype, only: dbl,i4b
use calc_func,only: an_r, an_phi, an_theta
  implicit none
  integer,intent(in) :: m
  integer,intent(out) :: k_r,k_phi,k_theta
  real(kind=dbl) :: mult
mult=1.d0
  ! part in r
  if((m.eq.1).or.(m.eq.9).or.(m.eq.17).or.(m.eq.25).or.(m.eq.3).or.(m.eq.11).or.(m.eq.19).or.(m.eq.27).or.&
       &(m.eq.4).or.(m.eq.12).or.(m.eq.20).or.(m.eq.28).or.(m.eq.7).or.(m.eq.15).or.(m.eq.23).or.(m.eq.31)) then
     k_r=1

  else if((m.eq.5).or.(m.eq.13).or.(m.eq.21).or.(m.eq.29).or.(m.eq.6).or.(m.eq.14).or.(m.eq.22).or.(m.eq.30)&
       &.or.(m.eq.2).or.(m.eq.10).or.(m.eq.18).or.(m.eq.26).or.(m.eq.8).or.(m.eq.16).or.(m.eq.24).or.(m.eq.32)) then
     k_r=2
mult=mult*an_r
  else if((m.eq.33).or.(m.eq.41).or.(m.eq.49).or.(m.eq.57).or.(m.eq.35).or.(m.eq.43).or.(m.eq.51).or.(m.eq.59)&
       &.or.(m.eq.36).or.(m.eq.44).or.(m.eq.52).or.(m.eq.60).or.(m.eq.39).or.(m.eq.47).or.(m.eq.55).or.(m.eq.63))then
     k_r=3

  else if((m.eq.34).or.(m.eq.42).or.(m.eq.50).or.(m.eq.58).or.(m.eq.37).or.(m.eq.45).or.(m.eq.53).or.(m.eq.61)&
       &.or.(m.eq.38).or.(m.eq.46).or.(m.eq.54).or.(m.eq.62).or.(m.eq.40).or.(m.eq.48).or.(m.eq.56).or.(m.eq.64))then
     k_r=4
mult=mult*an_r
  end if
  !part in phi
  if((m.eq.1).or.(m.eq.9).or.(m.eq.33).or.(m.eq.41).or.(m.eq.2).or.(m.eq.10).or.(m.eq.34).or.(m.eq.42).or.(m.eq.4)&
       &.or.(m.eq.12).or.(m.eq.36).or.(m.eq.44).or.(m.eq.6).or.(m.eq.14).or.(m.eq.38).or.(m.eq.46)) then
     k_phi=1

  else if((m.eq.3).or.(m.eq.11).or.(m.eq.35).or.(m.eq.43).or.(m.eq.5).or.(m.eq.13).or.(m.eq.37).or.(m.eq.45).or.&
       &(m.eq.7).or.(m.eq.15).or.(m.eq.39).or.(m.eq.47).or.(m.eq.8).or.(m.eq.16).or.(m.eq.40).or.(m.eq.48)) then
     k_phi=2
mult=mult*an_phi
  else if((m.eq.17).or.(m.eq.25).or.(m.eq.49).or.(m.eq.57).or.(m.eq.18).or.(m.eq.26).or.(m.eq.50).or.(m.eq.58).or.&
       &(m.eq.20).or.(m.eq.28).or.(m.eq.52).or.(m.eq.60).or.(m.eq.22).or.(m.eq.30).or.(m.eq.54).or.(m.eq.62))then
     k_phi=3

  else if((m.eq.19).or.(m.eq.27).or.(m.eq.51).or.(m.eq.59).or.(m.eq.21).or.(m.eq.29).or.(m.eq.53).or.(m.eq.61).or.&
       &(m.eq.23).or.(m.eq.31).or.(m.eq.55).or.(m.eq.63).or.(m.eq.24).or.(m.eq.32).or.(m.eq.56).or.(m.eq.64))then
     k_phi=4
mult=mult*an_phi
  end if
  !part in theta 

  if((m.eq.1).or.(m.eq.17).or.(m.eq.33).or.(m.eq.49).or.(m.eq.2).or.(m.eq.18).or.(m.eq.34).or.(m.eq.50).or.(m.eq.3)&
       &.or.(m.eq.19).or.(m.eq.35).or.(m.eq.51).or.(m.eq.5).or.(m.eq.21).or.(m.eq.37).or.(m.eq.53)) then
     k_theta=1

  else if((m.eq.4).or.(m.eq.20).or.(m.eq.36).or.(m.eq.52).or.(m.eq.6).or.(m.eq.22).or.(m.eq.38).or.(m.eq.54).or.(m.eq.7)&
       &.or.(m.eq.23).or.(m.eq.39).or.(m.eq.55).or.(m.eq.8).or.(m.eq.24).or.(m.eq.40).or.(m.eq.56)) then
     k_theta=2
mult=mult*an_theta
  else if((m.eq.9).or.(m.eq.25).or.(m.eq.41).or.(m.eq.57).or.(m.eq.10).or.(m.eq.26).or.(m.eq.42).or.(m.eq.58).or.(m.eq.11)&
       &.or.(m.eq.27).or.(m.eq.43).or.(m.eq.59).or.(m.eq.13).or.(m.eq.29).or.(m.eq.45).or.(m.eq.61))then
     k_theta=3

  else if((m.eq.12).or.(m.eq.28).or.(m.eq.44).or.(m.eq.60).or.(m.eq.14).or.(m.eq.30).or.(m.eq.46).or.(m.eq.62).or.(m.eq.15)&
       &.or.(m.eq.31).or.(m.eq.47).or.(m.eq.63).or.(m.eq.16).or.(m.eq.32).or.(m.eq.48).or.(m.eq.64))then
     k_theta=4
mult=mult*an_theta
  end if
end subroutine param_integr



!---------------------------
! Old integrand function
!---------------------------

real(kind=8) function func(x,y,zeta)
  use Calc_func
  use control
  use open_information
  use nrtype,only : dp,order,i4b
  use over, only : overlap_c
  use bspline
  implicit none
  integer(kind=i4b),parameter :: kxord=order,kyord=order,kzord=order
  real(kind=dp) basis,xx,yy,zz,real 
  real(kind=dp) :: x,y,zeta
  real(kind=dp) :: cos_phi,sin_phi,radius,rad,sin_theta,x2,x4,val 
  !  call c_func
  if (choice_matrix.eq.'gamma') then
     if (choice_partition.eq.'closed-closed') then 


        !----------------------------------------------------------------------
        !Using Green's theorem

        func=(radius(x)*der_r(x,k_ri)*der_r(x,k_rj)*bas_phi(y,k_phii)*bas_phi(y,k_phij)*bas_theta(zeta,k_thetai)&
             &*bas_theta(zeta,k_thetaj)/(an_r*an_r)+&
             &der_phi(y,k_phii)*der_phi(y,k_phij)*bas_r(x,k_ri)*bas_r(x,k_rj)*bas_theta(zeta,k_thetai)*bas_theta(zeta,k_thetaj)/&
             &(an_phi*an_phi*radius(x)*sin(angle_theta(zeta))*sin(angle_theta(zeta)))+&
             &der_theta(zeta,k_thetai)*der_theta(zeta,k_thetaj)*bas_r(x,k_ri)*bas_r(x,k_rj)*bas_phi(y,k_phii)*bas_phi(y,k_phij)/&
             &(an_theta*an_theta*radius(x)))*sin(angle_theta(zeta))*radius(x)*an_r*an_theta*an_phi
        !func =1
        !----------------------------------------------------------------------

     else if (choice_partition.eq.'open-closed') then
        func=((radius(x)*der_r(x,k_ri)*der_r(x,k_rj)*bas_phi(y,k_phii)*bas_phi(y,k_phij)*bas_theta(zeta,k_thetai)*&
             &bas_theta(zeta,k_thetaj)/(an_r*an_r)+&
             &der_phi(y,k_phii)*der_phi(y,k_phij)*bas_r(x,k_ri)*bas_r(x,k_rj)*bas_theta(zeta,k_thetai)*bas_theta(zeta,k_thetaj)/&
             &(an_phi*an_phi*radius(x)*sin(angle_theta(zeta))*sin(angle_theta(zeta)))+&
             &der_theta(zeta,k_thetai)*der_theta(zeta,k_thetaj)*bas_r(x,k_ri)*bas_r(x,k_rj)*bas_phi(y,k_phii)*bas_phi(y,k_phij)/&
             &(an_theta*an_theta*radius(x)))*sin(angle_theta(zeta))*radius(x))*an_r*an_theta*an_phi*&
             &open_vect(overlap_c(PatchVertex(element,beta),1),m)


        !func=1

     else if (choice_partition.eq.'open-open') then
        sin_theta=sin(angle_theta(zeta))
        rad=radius(x)
        x2=x*x
        x4=x2*x2
        func=((rad*der_r(x,k_ri)*der_r(x,k_rj)*bas_phi(y,k_phii)*bas_phi(y,k_phij)*bas_theta(zeta,k_thetai)*&
             &bas_theta(zeta,k_thetaj)/(an_r*an_r)+&
             &bas_r(x,k_ri)*bas_r(x,k_rj)/rad*(der_phi(y,k_phii)*der_phi(y,k_phij)*bas_theta(zeta,k_thetai)*bas_theta(zeta,k_thetaj)/&
             &(an_phi*an_phi*sin_theta*sin_theta)+&
             &der_theta(zeta,k_thetai)*der_theta(zeta,k_thetaj)*bas_phi(y,k_phii)*bas_phi(y,k_phij)/(an_theta*an_theta)))*sin_theta*rad)&
             &*an_r*an_theta*an_phi*open_vect(overlap_c(PatchVertex(rad_elem,beta),1),mm)*&
             &open_vect(overlap_c(PatchVertex(rad_elem,alpha),1),m)


        !func =1

        !        write(335,*)an_r,an_phi,der_ri(x,alpha),bas_ri(x,alpha)

     end if
  else if (choice_matrix.eq.'overlap') then
     if (choice_partition.eq.'closed-closed') then 
        func=(bas_r(x,k_ri)*bas_r(x,k_rj)*bas_phi(y,k_phii)*bas_phi(y,k_phij)*bas_theta(zeta,k_thetai)*bas_theta(zeta,k_thetaj))&
             &*radius(x)*sin(angle_theta(zeta))*radius(x)*an_r*an_theta*an_phi
        !func =1

     else if (choice_partition.eq.'open-closed') then

        func=((bas_r(x,k_ri)*bas_r(x,k_rj)*bas_phi(y,k_phii)*bas_phi(y,k_phij)*bas_theta(zeta,k_thetai)*bas_theta(zeta,k_thetaj))&
             &*radius(x)*sin(angle_theta(zeta))*radius(x))*an_r*an_theta*an_phi*open_vect(overlap_c(PatchVertex(element,beta),1),m)
        !func=1

     else if (choice_partition.eq.'open-open') then

        func=((bas_r(x,k_ri)*bas_r(x,k_rj)*bas_phi(y,k_phii)*bas_phi(y,k_phij)*bas_theta(zeta,k_thetai)*bas_theta(zeta,k_thetaj))&
             &*radius(x)*sin(angle_theta(zeta))*radius(x))*an_r*an_theta*an_phi
        !*open_vect(overlap_c(PatchVertex(rad_elem,beta),1),mm)*open_vect(overlap_c(PatchVertex(rad_elem,alpha),1),m)



        !func=1
     end if
  end if
end function func

function if_stack(r,theta,phi,epsilon_c,pos,angle,angle_phi)
  use nrtype, only : dbl,dpc,i4b,Pi
  logical :: if_stack
  real(kind=dbl),intent(in) :: r,theta,phi,epsilon_c,pos,angle,angle_phi
  !1- Normal atoms
  if ((  ((pos.ge.(r-epsilon_c)).and.(pos.le.(r+epsilon_c)))&
       &.and.((angle.ge.(theta-epsilon_c)).and.(angle.le.(theta+epsilon_c)))&
       &.and.((angle_phi.ge.(phi-epsilon_c)).and.(angle_phi.le.(phi+epsilon_c))) )&
       !2- atoms at center of grid
       !&.or.(  ((pos.ge.(r-epsilon_c)).and.(pos.le.r+epsilon_c))&
       !&.and.((angle.ge.(2*Pi-epsilon_c)).and.(angle.le.(2*Pi+epsilon_c)))  )&

       !3- atoms on z-axis 
     &.or.(  ((pos.ge.(r-epsilon_c)).and.(pos.le.(r+epsilon_c)))&
          &.and.((angle.ge.(theta-epsilon_c)).and.(angle.le.(theta+epsilon_c)))&
          &.and.((angle_phi.ge.(4*Pi-epsilon_c)).and.(angle_phi.le.(4*Pi+epsilon_c)))  )) then
     if_stack=.TRUE.
  else
     if_stack=.FALSE.
  end if
end function if_stack
