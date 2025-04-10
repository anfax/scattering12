!!!*************************************************************
! 文件/File: common_blocks.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: common_blocks.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

module precision
  ! set working precision
  ! this is used by all programs to be compatible with those 
  ! that remain in fortran77
  implicit none
  integer, parameter :: wp = kind(1.0D0)
end module precision
!----------------------------------------------------------------------------

module shared
  use precision
  implicit none 
  ! conversion constants
  real(kind=wp), parameter :: amu_to_au   = 0.182288853e+4_wp
  real(kind=wp), parameter :: cm_to_au    = 4.55633538e-6_wp
  real(kind=wp), parameter :: au_to_angst = 0.52917726_wp
  real(kind=wp), parameter :: pi          = 3.14159265358979_wp

  ! energy
  real(kind=wp) :: total_energy, reduced_e, two_mu, internal_energy

  ! propagator parameters
  real(kind=wp) :: rbegin, rswitch, rend 
  real(kind=wp) :: step_size, half_step_size, irxset          
  logical       :: use_variable_phase
  integer       :: nsteps     
  
  ! various
  logical       :: homonuclear          
  integer       :: rot_state_symm       
  integer       :: debug               
  logical       :: care
 
contains

  ! this is just a function to determine the parity of a positive integer
  ! on exit parity=.true. means an even number, parity=.false. means odd
  function parity(integer_number)
    integer :: integer_number, iparity
    logical :: parity
    if (integer_number<0) then
       write(*,*)'error in parity function: argument must be positive'
       write(*,*)'we were called with ', integer_number,' as argument'
       stop
    endif
    iparity=integer_number-2*(integer_number/2)
    if (iparity==1) parity=.false.
    if (iparity==0) parity=.true.
  end function parity

end module shared
!-----------------------------------------------------------------------------

module potential
  use precision
  ! save attribute for allocatable arrays is not needed since these
  ! are allocated in main
  implicit none
  integer :: num_lambda_terms
  integer :: num_pot_points 
  integer :: points_to_discard 
  integer :: lambda_to_zero
  ! Here there are memory segments for potential expansion.
  ! -v_lambda             contains vlambdas 
  ! -spline_coeff         contains the spline vectors for vlambdas
  ! -r                    is the R points vector
  ! -extrapolation_coeff  contains the extrapolation/fitting parameters
  !                         (a, b, c, d), its main dimension is 4.  
  real(kind=wp), allocatable, dimension(:,:,:,:) :: v_lambda, spline_coeff
  real(kind=wp), allocatable, dimension(:)      :: r
  real(kind=wp), allocatable, dimension(:,:,:,:) :: extrapolation_coeff
  ! this is to be passed to vpot
  integer :: vib_levels
end module potential

!------------------------------------------------------------------------------

module potential_scratch
  use precision
  ! these are constant size scratch vectors used in vpot
  !  allocate them in main together with potential module 
  !  allocatables (module just above)
  real(kind=wp), save, allocatable, dimension(:,:,:) :: vtt
  real(kind=wp), save, allocatable, dimension(:) :: v_tmp, v_tmp2 
end module potential_scratch

!------------------------------------------------------------------------------

module molecular_parameters
  use precision
  implicit none
  ! Here are molecular parameters.
  ! e_level is the energy of the quantum level
  ! rot_constant is the rotational constant associated to the vibrational level
  ! vibr_energy is the vibrational energy
  ! n_qn and j_qn are the quantum numbers associated with the
  !   quantum levels
  real(kind=wp),allocatable,dimension(:) :: e_level
  real(kind=wp),allocatable,dimension(:) :: rot_constant 
  real(kind=wp),allocatable,dimension(:) :: vibr_energy
  integer,allocatable,dimension(:) :: max_j_for_n, n_qn, j_qn, r_qn, s_qn 
  real(kind=wp) :: spin_spin  !spin-spin coupling constant in 1/cm
  real(kind=wp) :: spin_rot   !spin-orbit coupling constant in 1/cm  
end module molecular_parameters
!------------------------------------------------------------------------------

module xs
  use precision
  ! Here is the cross section that depends on the initial quantum levels
  !   and therefore is a rank 6 array.
  real(kind=wp), allocatable, dimension(:,:,:,:,:,:) :: total_xs2
  ! Partial contains partial cross sections
  !   they are initalized every jtot cycle. 
  real(kind=wp), allocatable, dimension(:,:,:,:,:,:) :: partial_xs2
  ! For each partial waves there are different symmetries
  !   just a variable to accumulate them
  real(kind=wp), allocatable, dimension(:,:,:,:,:,:) :: symmetry_xs2
  real(kind=wp), allocatable, dimension(:,:,:,:,:,:,:) :: sigma_J2
  real(kind=wp), allocatable, dimension(:) :: sigma_JJ2
  real(kind=wp), allocatable, dimension(:,:) :: xxtmp2 
  real(kind=wp), allocatable, dimension(:) :: sigma_el2
  real(kind=wp), allocatable, dimension(:) :: sigma_in2

end module xs
!------------------------------------------------------------------------------

module long_range_mod
  use precision
  ! These are the vectors containing the exponent for the long range
  ! behaviour of the v_lambdas.
  integer, allocatable, dimension(:) :: n_exp, m_exp
  character(len=2) :: select_propagator
end module long_range_mod
!-------------------------------------------------------------------------------

module var_ph_module
  use precision
  ! these variables are for var-ph propagator
  ! flag to select modified VP algorithm
  logical :: beszero, onlymod, realk
  ! some namelist variables
  real(kind=wp) :: eps, soglia, sigch
end module var_ph_module


module percival_seaton
  use precision
  ! potential_coefficient is (num_lambda_terms,local_channels,local_channels)
  ! this is allocated in close_coupling()
  real(kind=wp), allocatable, dimension(:,:,:) :: pot_coeff  
end module percival_seaton

!-------------------------------------------------------------------------------
module channels
  use precision
  ! local quantum numbers and energies
  ! these are allocated and used inside close_coupling()
  real(kind=wp), allocatable, dimension(:) :: local_energy 
  real(kind=wp), allocatable, dimension(:) :: wave_vector  
  integer, allocatable, dimension(:) :: j_local, n_local, &
                                        l_local, r_local, s_local 
end module channels








