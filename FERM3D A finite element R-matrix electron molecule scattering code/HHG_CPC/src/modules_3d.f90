!!!*************************************************************
! 文件/File: modules_3d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: modules_3d.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


MODULE nrtype
  implicit none
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I4 = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  integer,parameter :: dbl = kind(1.0d0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
  INTEGER, PARAMETER :: LGT = KIND(.true.)
  REAL(KIND=DBL), PARAMETER :: PI=3.141592653589793238462643383279502884197_dp
  REAL(KIND=DBL), PARAMETER ::PIO2=1.57079632679489661923132169163975144209858_dp
  REAL(KIND=DBL), PARAMETER ::TWOPI=6.283185307179586476925286766559005768394_dp
  REAL(KIND=DBL), PARAMETER ::SQRT2=1.41421356237309504880168872420969807856967_dp
  REAL(KIND=DBL), PARAMETER ::EULER=0.5772156649015328606065120900824024310422_dp
  REAL(KIND=DBL), PARAMETER ::PI_D=3.141592653589793238462643383279502884197_dp
  REAL(KIND=DBL), PARAMETER ::PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(KIND=DBL), PARAMETER ::TWOPI_D=6.283185307179586476925286766559005768394_dp
  real(kind=dbl), parameter ::nep_e=2.718281828459045090795598298427648842335_dp
  TYPE sprs2_sp
     INTEGER(I4B) :: n,len
     REAL(KIND=DBL), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp
  TYPE sprs2_dp
     INTEGER(I4B) :: n,len
     REAL(KIND=DBL), DIMENSION(:), POINTER :: val
     INTEGER(I4B), DIMENSION(:), POINTER :: irow
     INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp
  integer(kind=i4b),parameter :: order=5
END MODULE nrtype


module Brep
  use nrtype,only : dbl,i4b
  real(kind=dbl),dimension(:),allocatable :: Delel_r,Delel_phi,Delel_theta
  integer(kind=i4b),dimension(:,:),allocatable :: nodes,PatchVertex
  real(kind=dbl),dimension(:),allocatable :: r_in,r_out,phi_in,phi_out,theta_in,theta_out,&
       &r_glo,phi_glo,theta_glo,r_glo_avg,phi_glo_avg,theta_glo_avg
  real(kind=dbl),dimension(:),allocatable :: U,coord_x,coord_y,coord_z
  integer(kind=i4b),dimension(:),allocatable :: nr_z,nphi_z,ntheta_z,zone_flag,nodes_closed_index
  integer(kind=i4b) :: n,p,p_Patch,z,phi_sectors,plotpoints
  real(kind=dbl),parameter :: epsilon_custom=1.0d-4
  real(kind=dbl) ::R0,shift
  character*10,dimension(:),allocatable :: boundary,edge
  character*10,dimension(:,:),allocatable :: full_boundary

  !-------------------------------------------------------------------------
  ! This module provides the grid arrays to the other subroutines and in 
  ! general info about the box 
  ! Variables:
  ! Delta_r,Delta_phi = step in r,phi, within the sectors (depends on zone)
  ! nodes = physical nodes indices 
  ! PatchVertex = index of the functions (stored as 8 fictitious nodes for each
  !   physical node), bridges local and global identification of the nodes
  ! r_in,r_out,phi_in,phi_out = first and last r,phi for each zone
  ! r_glo,phi_glo,theta_glo = r,phi,theta for the physical nodes
  ! phi_glo_in = 
  ! r_glo_avg,phi_glo_avg = average r,phi within each sector
  ! U = potential
  ! nr_z,nphi_z = number of sectors within a zone
  ! zone_flag = number of zone in which a function is 
  ! n = number of sectors
  ! p = number of physical nodes
  ! p_Patch = number of functions (fictitious nodes)
  ! shift = shift for z-axis to put atom with most of charge on origin, to avoid
  ! dipole field in eteronuclear diatomics
  ! z =  number of zones
  ! R0 = dimension of the box
  ! boundary = position of each physical node (bulk, up, down,left,right)
  ! edge = position of each function
  ! full_boundary = same as edge, but with two indices (element,#)
  ! zone_flag = identifies the zone number for scaling purposes in the integrals
  !-------------------------------------------------------------------------
end module Brep

module control
  use nrtype, only : i4b
  integer(kind=i4b),parameter :: problem=2
  integer(kind=i4b) :: partition,choice_in_func,output_file
  character*15 :: choice_matrix,choice_open,choice_func,choice_partition,control_integr  
  character*10,parameter :: restore_gamma='no',restore_overlap='no',iter_cont='yes',&
       & option_solver='Pardiso',option_Coulomb='Seaton',option_diagonalizer='Arpack'
  !,molecule='neutral'
  character*10 :: restore='no',option_wf='nowf',&
       & calculation_type='scatter',exchange='slater',molecule='neutral'
  character(len=15) :: option_gauge='acceleration'
  logical :: cont=.TRUE.
  logical :: DFT=.FALSE.
  !-------------------------------------------------------------------------
  ! Module that bears the control variables
  ! problem = obsolete
  ! partition = appears in the external do loop in mat
  ! choice_matrix = choice between 'gamma' and 'overlap'
  ! choice_open = choice between 'open' and 'overlap' in overlap
  ! choice_func = choice between 'pot','func_en'
  ! choice_partition = choice between 'closed-closed','open-closed','open-open'
  ! control_integr = sing or outside (distinguishes functions that have value !=
  !   0 on the singularity for integration purposes)
  ! calculation_type = scatter or bound
  ! exchange = slater or hara, in this version hara exchange is linearized in
  !   energy so it is called with the "slater" option and slater exchange cannot be
  !   called
  ! molecule = ion or neutral
  ! iter_cont = "yes" chooses an iterative algorithm to solve linear system, "no"
  !   selects SuperLU direct methods 
  !-------------------------------------------------------------------------
end module control


module Open_information
  use nrtype,only : dbl,i4b
  integer(kind=i4b),dimension(:,:),allocatable :: open_info,long_open
  integer(kind=i4b),dimension(:),allocatable :: wf_ind,rad_elem_v,alpha_v
  real(kind=dbl),dimension(:,:),allocatable :: open_vect
  integer(kind=i4b) :: max_open,max_index,alpha,beta,counter,rad_elem,lmax,long_max,vect_max,max_open_total
  character*10,dimension(:,:),allocatable :: bounda_boundary,bounda2_boundary
  character*10,dimension(:),allocatable :: bionda2_boundary,bionda_boundary
  real(kind=dbl) :: kappa,E,Energy_of_cycle
  integer(kind=i4b) :: max_l=8
  save max_open
  save max_index
  !-------------------------------------------------------------------------
  ! Brings info about harmonics on the surface
  ! open_info = function,element and # within the sector for each open function
  ! open_vect = coefficient of the open function
  ! max_open = max number of open channels 
  ! max_index = max number of open functions
  ! alpha = # of function within the sector for the open partitions
  ! beta = same as alpha
  ! counter,counterx = counters
  ! rad_elem = # of element for the open partitions 
  ! bounda_boundary / bounda2_boundary = location of the function (for circular boundary): 'up_fun'
  !                   or 'none',2 indices
  ! bionda_boundary / bionda2_boundary = same as before, 1 index
  !-------------------------------------------------------------------------
end module Open_information

module Matrices
  use nrtype,only : dbl,i4b
  real(kind=dbl),dimension(:,:),allocatable :: gamma,area_overlap,open_open,open_closed,area_overlap_open_open,area_overlap_open_closed,d,open_mat,overlap,gamma_temp
  integer(kind=i4b),dimension(:,:),allocatable :: gamma_c,gamma_c_temp
  integer(kind=i4b),dimension(:),allocatable :: area_overlap_c
  integer(kind=i4b) :: nnz
  !-------------------------------------------------------------------------
  ! gamma = del*del matrix for closed-closed partition
  ! area_overlap = psi*psi for closed-closed partition
  ! open_open = del*del matrix for open-open partition
  ! open_closed = del*del for open-closed partition
  ! area_overlap_open_open = psi*psi matrix for open-open partition
  ! area_overlap_open_closed = psi*psi for open-closed partition
  ! gamma_c = stores index pattern for rows in sparse format (closed-closed
  !    partition)
  ! area_overlap_c = stores the index number f a function just among the
  !    closed ones
  !-------------------------------------------------------------------------
end module Matrices


module Calc_func
  use nrtype,only : dbl,i4b
  use Brep
  use Open_information
  use control
  implicit none
  integer(kind=i4b) :: element,i,j,index,indexx,lj,mj,fj,li,mi,fi,m,mm,op_index,op_indexx
  real(kind=dbl) :: an_r,an_phi,dn_r,an_theta
  integer(kind=i4b) :: k_ri,k_rj,k_thetai,k_thetaj,k_phii,k_phij
  real(kind=dbl),EXTERNAL ::bas_r,bas_phi,bas_ri,bas_rj,bas_phii,bas_phij,der_ri,&
       &der_rj,der_phii,der_phij,sec_derbas,der_r,der_phi,der_bas,&
       &sec_derbas_rj,sec_derbas_phij,bas_thetai,bas_thetaj,der_thetai,&
       &der_thetaj,sec_derbas_thetaj,angle_theta,bas_theta,der_theta,angle_phi
  real(kind=dbl),dimension(:,:), allocatable :: interp_V,interp_V_sing,interp_density,interp_density_sing
  real(kind=dbl),dimension(:,:), allocatable :: interp_grad,interp_grad_sing,interp_nabla,interp_nabla_sing
  real(kind=dbl),dimension(:,:), allocatable :: der_tab, bas_tab, der_tab_sing, bas_tab_sing
  integer(kind=i4b) :: num,ind_sing
  integer(kind=i4b),dimension(:),allocatable :: vect_sing
  integer(kind=i4b) :: count_func,count_pts
  real(kind=dbl),dimension(:,:),allocatable :: mem_tab
  real(kind=dbl),dimension(:),allocatable :: val_array, val_array_sing
  real(kind=dbl) :: Energy_Hara,polar,polar_2,r_cutoff,Energy_last_bound
  !-------------------------------------------------------------------------
  ! k_ri,k_rj,k_thetai,k_thetaj,k_phii,k_phij = indices for the polynomials in r,
  !   phi, theta depending on the local index of the function i or j inside the
  !   element  
  ! der_tab, bas_tab, der_tab_sing, bas_tab_sing = arrays that store the values
  !   of the polynomials and derivatives for each integration point in an element
  ! element,i,j,index,indexx,m,mm = indices
  ! an_r,an_phi = scale factor for the sector (width of the sector) for
  !               integration
  ! dn_r = radius of each node in the sector
  ! Energy_last_bound = Energy of the last bound state of the molecule (input to Hara potential)
  ! r_cutoff = cutoff radius for polarization 
  ! mem_tab = array that stores the value of the function products for the
  !   integrals
  !-------------------------------------------------------------------------
contains
  subroutine c_func
    use Brep
    use control
    if (choice_partition.ne.'open-open') then
       an_r=Delel_r((element))
       an_phi=Delel_phi((element))
       an_theta=Delel_theta((element))
    else 
       an_r=Delel_r((rad_elem))
       an_phi=Delel_phi((rad_elem))
       an_theta=Delel_theta((rad_elem))
    end if
  end subroutine c_func
end module Calc_func




module bsp_info
  ! Information for the bsplibs for parallelization of the algorithm
  use nrtype, only : i4b
  integer(kind=i4b) :: bsppid, bspnprocs,nE
  integer(kind=i4b), parameter :: bspint=4,bspreal=4,bspdouble=8
  !!external bsppid, bspnprocs
end module bsp_info
