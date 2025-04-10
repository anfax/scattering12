!!!*************************************************************
! 文件/File: main.cpp.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: main.cpp.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

! Program written by D. Lopez-Duran, E. Bodo, and F. A. Gianturco.
! April 2008, version 0408.
! 
! This program solves the close coupling equations for atom-vibrating
!   diatom scattering. This version of main() uses mpi calls to parallelize
!   the jtot loop and #ifdef preprocessor directives to include the non-mpi
!   serial version.
! 
! To understand the physical approach of the problem see the next
!  references:
!  1) Russel T. Pack, J. Chem. Phys. 60, 633-639 (1974) (in this
!     article spin equal to 0, in this program spin equal or not to 0).
!     This article will be referred in this code as "Pack".
!  2) A. M. Arthurs and A. Dalgarno, Proceedings of the Royal Society of
!     London. Series A, Mathematical and Physical Sciences, Vol. 256,
!     No. 1287 (Jul. 19, 1960), pp. 540-551.
!
! The main program reads from different namelists the input parameters
!  and from an external unit the vibrational coupling elements as a 
!  function of r for each gamma, n and n'.
!
! The input parameters are commented in the code at the read/namelist/
!  statements.
!
! This program is valid in both homonuclear and heteronuclear diatoms, and
!  in any state of spin (spin:= singlet, doublet, or triplet).
! The orbital angular momentum of the diatom is 0, i.e., we are in a sigma
!  state of the diatom.
!
! Vectors and quantum numbers:
!  N:= pure rotational angular momentum of the diatom
!  s:= diatom spin
!  j:= N+s, total angular momentum of the diatom
!  (not of the whole system!)
!  j, N, and s satisfy the triangular inequality.
!  v:= vibrational quantum number
!
! Arrays used in the code to simulate this vectors:
!  n -> v, vibrational quantum number
!  j -> j, total angular momentum quantum number of the diatom
!  (not of the whole system!)
!  r -> N, pure rotational angular momentum quantum number
!
! The structure of the program is quite simple:
!
!  main:  1.  root process examines input,
!             (send all input datas to other process)
!         2.  sets correct parity cases for homo. molec.,
!         3.  counts and orders asymptotic channels,
!         4.  root peocess reads v_lambda,
!             (send potential to other process)
!         5.  does a cycle over energies, 
!         5.5 estabish which jtot values will be done on each mpi process
!             depending on the number of processes chosen. Send the 
!             choosen values
!             to all other processes  
!         5.6 every process starts its loop 
!         6.  inside the loop, each process:
!             calls  subroutines to handle cc equations
!         7.  accumulates total-xs,
!         7.5 root process collect al data through mpi_reduce()
!         8.  writes output.
!
program scattering

  !+---------------------+
  !| use various modules |
  !+---------------------+
  ! ->Set working precision to old-F77 double precision
  use precision

  ! ->Conversion factors and constants  
  ! ->Total energy(E), 2*mu*E, 2*mu
  ! ->Propagator parameters
  ! ->Debugging level of the program
  ! ->Variables for the long_range propagator choice
  use shared

  use molecular_parameters

  ! ->Cross sections  
  use xs

  use potential
  
  use potential_scratch

  use long_range_mod
  
  ! ->Accuracy parameters for var_ph propagator
  use var_ph_module
  !--------------------------------------------

  implicit none

#ifdef PARALLEL
  include "mpif.h"
#endif

  !+-----------------+
  !| interface block |
  !+-----------------+
  ! this is needed for the xs() assumed shape array in cc.f   
  interface
     subroutine close_coupling(jtot, step_par, asympt_channels, &
          open_channels, local_channels, symmetry, xs2,&
          myrank, nproc, first_time, multiplicity)
       real(kind=kind(1.0D0)),dimension(:,:,:,:,:,:),intent(OUT) :: xs2
       integer, intent(IN)  :: asympt_channels
       integer, intent(INOUT)  :: local_channels 
       integer, intent(IN)  :: jtot
       integer, intent(IN)  :: step_par
       integer, intent(IN)  :: symmetry
       integer, intent(OUT) :: open_channels
       integer, intent(in)  :: myrank, nproc
       logical, intent(IN)  :: first_time
       integer, intent(IN)  :: multiplicity
     end subroutine close_coupling
  end interface
  !------------------------------------------

  !+-----------------------------------------------------+
  ! this is for compatibility with the old xf3j and xf6j.|
  ! just some factorial for those subroutines            |
  ! ## check if we can change ##                         |
  real(kind=wp) :: fact                                  !
  common/fct/fact(10000)                                 !
  !+-----------------------------------------------------+ 
  
  !+-------------------+
  !| declaration block |
  !+-------------------+
#ifdef PARALLEL
  !+-------------- parallel enviroment --------------+
  integer :: myrank     ! rank of the process 
  integer :: nproc      ! number of mpi processes
  integer :: ierr_mpi   ! error code from mpi calls
  integer :: dim1, dim2
#endif

#ifndef PARALLEL
  ! These are needed to simulate mpi values inside cc(), vfitt(),
  !   vpot() and daprop() calls for correct output handling
  integer, parameter :: myrank=0, nproc=1
  integer, dimension(1) :: jtotm, jtotx
#endif

  ! fortran units in which we put every write(*,*)
  ! statement that is controlled by a debug>=1 condition
  ! for nproc==1 (both mpi or not) unit_mpi becomes std-output (unit=6)
  integer :: unit_mpi 

#ifdef PARALLEL
  ! various number of partial waves
  integer :: n_pw, pw_per_proc, pw_remaining, jtot_temp

  ! vectors containing the values of jtot to be done on each process
  integer, allocatable, dimension(:) :: jtotm, jtotx, jtot_per_proc

  ! cumulative rank0 cross section to build an mpi_reduce() call
  ! (this is really allocated only by root process)
  real(kind=wp), allocatable, dimension(:,:,:,:,:,:) :: real_total2
#endif  
  !+------------------------------------------------------------------+

  ! these are needed by the relaxed=.true. case
  real(kind=wp) :: e_new
  integer       :: n_new, j_new
  integer       :: r_new, s_new 

  !These are for variables that are in namelists, 
  ! but not declared in commons
  real(kind=wp) :: reduced_mass 
  real(kind=wp) :: first_energy,delta_e  
  integer :: num_vib_levels, max_j_value
  integer :: initial_vib_level, init_vibrat_level     
  integer :: num_of_energies, steps=20 , init_N_state
  integer :: potential_unit ! unit from which main() reads the v_lambda  
  real(kind=wp), dimension(100) :: b, e_vib, xjmax 
  integer      , dimension(100) :: jmax, nlr, mlr
  real(kind=wp) :: xinit_j_state, xjout, xmax_j_value
 
  !These are local variables
  real(kind=wp) :: e, e_to_compare, collision_energy
  real(kind=wp) :: rend_backup,rswitch_backup
  integer :: num_rot_levels
  integer :: num_N_levels  
  integer :: omega, omega_min, omega_max     
  integer :: initial_j, j_step ! j_parameters for homo molecules       
  integer :: initial_r 
  integer :: i,j,l,n,kr,ii,nn
  integer :: ir               
  integer :: en_index         
  integer :: count            
  logical :: strange_case=.false. ! only for parallel
  integer :: jtemp
  logical :: logarithmic_e, first_time, relaxed
  logical :: singlet, doublet, triplet
  integer :: multiplicity, ccount
  real(kind=wp) :: corr
  character (len=7) :: stsp

  ! Local channels label the |j,l,n,r> states 
  ! while asymptotic channels label the |j,n,r> states.
  integer :: num_ch         ! total number of asymptotic channels(n,j,r)
  integer :: num_local_ch   ! total number of local channels(n,j,l,r)
  integer :: open_ch        ! open asymptotic channels(n,j,r)

  real(kind=wp) :: v, xtmp2
  real(kind=wp) :: dummy   
  real(kind=wp) :: prec 
  integer :: ierr              ! err code from i/o operation
  integer :: j1,j2,n1,n2,r1,r2,s1,s2 
  integer :: choice           
  integer :: jtot, jtotmin,jtotmax, jtotstep
  real(kind=wp) :: xjtot, xjtotmin, xjtotmax, xjtotstep
  integer :: jreal, lmin, lmax            
  integer :: rreal  
  integer :: term1, term2, term3, term4, term5, term6, term7
  integer :: uno, dos, tres
  integer :: i_dummy          

  namelist /molecular_general/ num_vib_levels, xmax_j_value, &
       init_vibrat_level, xinit_j_state, init_N_state, xjout, & 
       reduced_mass, relaxed, singlet, doublet, triplet
  namelist /energy/ first_energy, delta_e, num_of_energies,logarithmic_e
  namelist /potential_general/ num_lambda_terms, num_pot_points, &
       points_to_discard, lambda_to_zero
  namelist /molecular_params/ b, xjmax, e_vib, homonuclear, &
       rot_state_symm, spin_spin, spin_rot
  namelist /program_behaviour/ debug, potential_unit
  namelist /propagator/ rbegin, rswitch, rend, steps, irxset, &
       use_variable_phase, eps, soglia, &
       sigch, onlymod, care
  namelist /partial_waves/ xjtotmin, xjtotmax, xjtotstep
  namelist /exponents/ nlr, mlr
  
  ! +--------------------------------+
  ! |    Code initialization done    |
  ! |    See file common_blocks.f90  |
  ! +--------------------------------+
  
#ifdef PARALLEL
  ! initialize the mpi subsystem
  call mpi_init(ierr_mpi)
  call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr_mpi)
  call mpi_comm_size(MPI_COMM_WORLD, nproc,  ierr_mpi)
#endif

  ! +---------------------------------------+
  ! |  rotational and vibrational channels  |
  ! |  namelist  molecular_general          |
  ! +---------------------------------------+
  ! -num_vib_levels is the total number of channels
  !    (indexing of vibrational channels start from 1 
  !     in all the sequent code)
  ! -init_vibrat_level is the initial vibrational state
  !    (0 means ground state)
  ! -xinit_j_state is the initial total rotational state (j=N+s)
  ! -init_N_state is the initial pure rotational (diatomic) state
  ! -max_j_value is the maximum j used in the vibrational folders: 
  !    (rotational channel are labelled from 0) 
  !    (jmax _must_ be larger than every element of max_j_for_n(:))
  ! -reduced_mass is the reduced mass in a.m.u
  !
  ! -possible states:
  ! spin=0   -> one only possible state
  ! spin=1/2 -> F1 and F2
  ! spin=1   -> F1, F2, and F3
  ! (see J. Brown and A. Carrington, Rotational Spectroscopy of
  !  Diatomic Molecules, pages 21 to 23, for instance)
  relaxed=.false.
  reduced_mass=1.0_wp
  num_vib_levels=5
  xmax_j_value=15.0_wp
  init_vibrat_level=0
  xinit_j_state=0.0_wp
  init_N_state=0
  singlet=.true.
  doublet=.false.
  triplet=.false.
  xjout=0.0_wp
  if (myrank == 0 ) then  !read input
     read(*,nml = molecular_general)
  endif

  initial_vib_level=init_vibrat_level+1

! consistency check on the vibrational levels 

   if(num_vib_levels<=init_vibrat_level)then
   write(*,*)'There are', num_vib_levels,'available vibrational levels' 
   write(*,*)'and you are requesting the initial vibrational &
              level number', init_vibrat_level+1,'.' 
   write(*,*)'(remember that the ground state is counted as 0)'
   write(*,*)''
   write(*,*)'It is wrong. Please check it.' 
   write(*,*)''
   write(*,*)'Program stop.'
   stop
  endif 

! consistency check on the angular momenta 

  if (xmax_j_value<xjout) then
   write(*,*)''
   write(*,898)'Max. j value=',xmax_j_value, &
               'jout=',xjout
   write(*,*)''
   write(*,*)'jout cannot be greater than the maximum j value. &
              Please check it.'
   write(*,*)''
   write(*,*)'Program stop.'
   stop
  endif  

  if(xinit_j_state>xmax_j_value)then
     write(*,*)''
     write(*,'(1x,a,1x,f4.1)')'max_j_value=',xmax_j_value
     write(*,'(1x,a,1x,f4.1)')'init_j_state=',xinit_j_state
     write(*,*)'But init_j_state cannot be bigger than max_j_value.'
     write(*,*)'Please check it.'
     write(*,*)''
     write(*,*)'Program stop.'
     stop
  endif

  ! we need also to know the multiplicity
  ! of the spin state of the diatom
  ! when calling close_coupling subroutine
  
   ccount=0
   if(singlet)then
   multiplicity=1
   stsp='singlet'
   corr=0.0_wp
   ccount=ccount+1
   endif
   if(doublet)then
   multiplicity=2
   stsp='doublet'
   corr=0.5_wp
   ccount=ccount+1
   endif
   if(triplet)then
   multiplicity=3
   stsp='triplet'
   corr=0.0_wp
   ccount=ccount+1
   endif

! consistency check on the spin state

   if(ccount==0.or.ccount>1)then
   write(*,*)''
   write(*,*)'Wrong spin state. Please check it.'
   write(*,*)''
   write(*,*)'Program stop.'
   stop
   endif

   if(singlet) then
     if(nint(xmax_j_value+0.4_wp)/=nint(xmax_j_value-0.4_wp)) then
     write(*,*)''
     write(*,899)'Max j value=',xmax_j_value
     write(*,*)'Non-physical choice of the maximum j value. &
                Please check it.'
     write(*,*)''
     write(*,*)'Program stop.'
     stop
     endif
     if(nint(abs(xjout+0.4_wp))/=nint(abs(xjout-0.4_wp))) then
     write(*,*)''
     write(*,899)'jout=',xjout
     write(*,*)'Non-physical choice of jout. Please check it.'
     write(*,*)''
     write(*,*)'Program stop.'
     stop
     endif
     if(abs(xinit_j_state-real(init_N_state,kind=wp))/=0) then
     write(*,*)''
     write(*,990)'initial j=',xinit_j_state, &
                 'initial N=',init_N_state
     write(*,*)'Non-physical choice of j and N. Please check it.'
     write(*,*)''
     write(*,*)'Program stop.'
     stop
     endif
   elseif(doublet) then
     if(nint(xmax_j_value+corr+0.4_wp)/=nint(xmax_j_value+corr-0.4_wp)) then
     write(*,*)''
     write(*,899)'Max j value=',xmax_j_value
     write(*,*)'Non-physical choice of the maximum j value. &
                Please check it.'     
     write(*,*)''
     write(*,*)'Program stop.'
     stop
     endif
     if(xjout<0.5) then
     write(*,*)''
     write(*,899)'jout=',xjout
     write(*,*)'Non-physical choice of jout. Please check it.'
     write(*,*)''
     write(*,*)'Program stop.'
     stop
     endif
     if(nint(abs(xjout+corr+0.4_wp))/=nint(abs(xjout+corr-0.4_wp))) then
     write(*,*)''
     write(*,899)'jout=',xjout
     write(*,*)'Non-physical choice of jout. Please check it.'
     write(*,*)''
     write(*,*)'Program stop.'
     stop
     endif
     if(abs(xinit_j_state-real(init_N_state,kind=wp))/=0.5_wp) then
     write(*,*)''
     write(*,990)'initial j=',xinit_j_state, &
                 'initial N=',init_N_state
     write(*,*)'Non-physical choice of j and N. Please check it.'
     write(*,*)''
     write(*,*)'Program stop.'
     stop
     endif
   elseif(triplet) then
     if(nint(xmax_j_value+0.4_wp)/=nint(xmax_j_value-0.4_wp)) then
     write(*,*)''
     write(*,899)'Max j value=',xmax_j_value
     write(*,*)'Non-physical choice of the maximum j value. &
                Please check it.'     
     write(*,*)''
     write(*,*)'Program stop.'
     stop
     endif
     if(nint(abs(xjout+0.4_wp))/=nint(abs(xjout-0.4_wp))) then
     write(*,*)''
     write(*,899)'jout=',xjout
     write(*,*)'Non-physical choice of jout. Please check it.'
     write(*,*)''
     write(*,*)'Program stop.'
     stop
     endif
     if(abs(real(xinit_j_state,kind=wp)-init_N_state)/=0.0_wp.and. &
        abs(xinit_j_state-real(init_N_state,kind=wp))/=1.0_wp.or. &
        xinit_j_state==0.0_wp.and.init_N_state==0) then
     write(*,*)''
     write(*,990)'initial j=',xinit_j_state, &
                 'initial N=',init_N_state
     write(*,*)'Non-physical choice of j and N. Please check it.'
     write(*,*)''
     write(*,*)'Program stop.'
     stop
   endif
   endif

   max_j_value=nint(xmax_j_value+corr)

#ifdef PARALLEL
  call mpi_bcast(num_vib_levels,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(max_j_value,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(initial_vib_level, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(init_vibrat_level, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(xinit_j_state, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(init_N_state, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(xjout, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(max_j_value, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(xmax_j_value, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(reduced_mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(relaxed, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(singlet, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(doublet, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(triplet, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(multiplicity, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(corr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
#endif

  ! +------------------------------+
  ! |  collision energies in cm-1  |
  ! |  namelist energy             |
  ! +------------------------------+
  ! # performing a cycle over energy value on a regular grid  #
  ! # the fact that this cycle will be the otermost one is    #
  ! # probably quite inefficient (we have to manage this).    #  
  ! -first_energy is collision energy in cm-1
  ! -delta_e is energy spacing in cm-1
  ! -number_of_energies is the number of energies
  !   e(i)=first_energy+deltae*(i-1), i=1,number_of_energies
  !
  ! logaritmic_e=.true. means that the collision energies are calculated
  ! on a logarithmic scale
  !   e(i)=(first_energy)*delta_e**(en_index-1)
  logarithmic_e=.false.
  first_energy=1000.0_wp
  delta_e=1.0_wp 
  num_of_energies =1
  if (myrank==0) then
     read(*,nml = energy)  
  endif
  
#ifdef PARALLEL
  call mpi_bcast(first_energy, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(delta_e,      1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(num_of_energies,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(logarithmic_e,   1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr_mpi)
#endif

  ! +-----------------------------+
  ! |  potential treatment input  |
  ! |  namelist potential_general |
  ! +-----------------------------+
  ! -num_pot_points is the number of R points in the couplings elements
  ! -num_lambda_terms is the number of terms in the potential expansion 
  ! -points_to_discard is used to start inverse power extrapolation
  !         in an inner r point and not at the end of the r array.
  !         the standard procedure begins extrapolation 
  !         at r(num_pot_points-1), 
  !         when points_to_discard.ne.0 the procedure will start at 
  !         r(num_pot_points-points_to_discard-1).
  ! -all the terms with lambda+1 > lambda_to_zero will be 
  !  setted to zero at long range
  !         es: if lambdas are 0,2,4,6 index goes through 1,2,3,4:
  !         if ilam= 2 then lambda 4 and 6 (index 3,4) will be 
  !         setted to zero 
  !         set ilam>=lamx if you don't want this option
  num_lambda_terms=0
  num_pot_points=0
  points_to_discard=0
  lambda_to_zero=0
  if (myrank==0) then
     read(*,nml = potential_general)
  endif

#ifdef PARALLEL
  call mpi_bcast(num_lambda_terms,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(num_pot_points,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(points_to_discard, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(lambda_to_zero,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
#endif

  ! +------------------------------------------+      
  ! |            molecular constants           |
  ! |       namelist molecular_params          |
  ! +------------------------------------------+      
  ! rot_constant(1:num_vib_levels) are the v-dependent 
  !  rotational const in cm-1
  !     (we read b(1:100)) 
  ! vibr_energy(1:num_vib_levels) are vibrational energy in cm-1:
  !  note: you can use either vibr_energy(1)=0.0d0 or vibr_energy(1)=zpe.
  !  (we read e_vib(1:100))   
  ! max_j_for_n(1,...,num_vib_levels) are the max. value of j for i
  !  each vib. manifold
  !     note: you are not requested to give a value for max_j_for_n
  !           just set them .lt.0 to use max_j_value+1 channels in each
  !           vibrational level.
  !     we read(jmax(1:100))   
  ! homonuclear is .true. for homonuclear diatoms.
  ! rot_state_symm is for homonuclear only:
  !     rot_state_symm set to 0 means we are using even j values
  !     rot_state_symm set to 1 means we are using odd  j values
  num_rot_levels=max_j_value+1  
  num_N_levels=(num_rot_levels)*3-2 
  allocate(max_j_for_n(num_vib_levels))
  allocate(rot_constant(num_vib_levels))
  allocate(vibr_energy(num_vib_levels))

  rot_constant=7.0_wp
  spin_spin=0.0_wp  
  spin_rot=0.0_wp   
  max_j_for_n=-1
  vibr_energy=1000.0_wp
  homonuclear=.false.
  rot_state_symm=0
  if (myrank==0) then
     read(*, nml = molecular_params)  
     do i=1,num_vib_levels
        jmax(i)=nint(xjmax(i)-corr)
        max_j_for_n(i)=jmax(i)
        rot_constant(i)=b(i)
        vibr_energy(i)=e_vib(i)
     enddo
  endif

! consistency check on the rotational state
  if (homonuclear) then
     if(rot_state_symm==0) then 
        if (.not.parity(init_N_state)) then
           write(*,*)''
           write(*,'(x,a,i3)')'initial_N_state = ',init_N_state
           write(*,*)''
           write(*,*)'But you have chosen *ODD* parity for & 
                      homonuclear molecule. Please check it.'
           write(*,*)''
           write(*,*)'Program stop.'
#ifdef PARALLEL
           write(*,*)'stopping mpi subsystem'
           call mpi_finalize(ierr_mpi)
#endif

           stop
        endif
     else
        if (parity(init_N_state)) then
           write(*,*)''
           write(*,'(x,a,i3)')'initial_N_state = ',init_N_state
           write(*,*)''
           write(*,*)'But you have chosen *EVEN* parity for & 
                      homonuclear molecule. Please check it.'
           write(*,*)''
           write(*,*)'Program stop.'

#ifdef PARALLEL
           write(*,*)'stopping mpi subsystem'
           call mpi_finalize(ierr_mpi)
#endif

           stop
        endif
     endif
  endif

! consistency check on the angular momenta
    if(singlet.or.triplet) then
     do i=1,num_vib_levels
       if(nint(xjmax(i)+0.4_wp)/=nint(xjmax(i)-0.4_wp)) then
         write(*,*)''
        write(*,897)'The vib level number',i,'has a jmax value of:',xjmax(i)
         write(*,*)'This is non-physical. Please check it.'
         write(*,*)''
         write(*,*)'Program stop.'
         stop
       endif
     enddo
    elseif(doublet) then
     do i=1,num_vib_levels
       if(nint(xjmax(i)+corr+0.4_wp)/=nint(xjmax(i)+corr-0.4_wp)) then
         write(*,*)''
        write(*,897)'The vib level number',i,'has a jmax value of:',xjmax(i)
         write(*,*)'This is non-physical. Please check it.'
         write(*,*)''
         write(*,*)'Program stop.'
         stop
       endif
     enddo
    endif

#ifdef PARALLEL
  dim1=num_vib_levels
  call mpi_bcast(max_j_for_n,  dim1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(rot_constant, dim1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(vibr_energy,  dim1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(homonuclear,    1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(rot_state_symm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)   
  call mpi_bcast(spin_rot,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(spin_spin,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
#endif

  ! +---------------------------------+
  ! |   namelist program_behaviour    | 
  ! +---------------------------------+
  ! debug is for debugging: set to 0 for normal calc., to 1 for 
  !     normal debug, to 2 to get information about details, 3 for 
  !     filling up the disk.
  !-------------------------------------------------------------------------+
  ! Debug is a little more difficult here with mpi calls and more than      |
  ! one process writing out. If you don't need any particular output you    |
  ! are supposed to use only debug=0 level. This level gives you the normal |
  ! output as in the serial version.                                        |
  ! If you need debug informations then take care that every write(*)       |
  ! statement controlled by a (debug>=1) is written on a fortran unit       |
  ! which number (#name#) is that of (20+myrank).                           |
  !-------------------------------------------------------------------------+
  ! potential_unit is the fortran unit number from which the program must
  !     read the v_lambda couplings.
  potential_unit=10
  debug=0
  if (myrank==0) then
     read(*,nml = program_behaviour)
  endif

#ifdef PARALLEL
  call mpi_bcast(debug,          1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(potential_unit, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
#endif
  
  ! +----------------------------+
  ! |   propagator parameters    |
  ! |   namelist propagator      |
  ! +----------------------------+
  ! rbegin is the initial point in integration; choose it somewhat 
  !     inside the classical forbidden region of the repulsive 
  !     potential. A value of 0.2 or 0.5*r(turning-point) is generally 
  !     resonable. (angstroms)
  ! rend is the matching point; be sure your potential is 
  !     10*-3 times the collision energy. If the centrifugal 
  !     barrier at long range is higher than open channels energy
  !     the program will stop, or, if you use the fitting method of
  !     potfitt rend will be increased. (angstroms)
  ! nsteps is the number of steps for the log-derivative propagator
  ! irxset setted to <0 means that the program will not try to evaluate
  !     a better rend for integration when some of the asymptotic channels 
  !     lie under the centrifugal barrier.irxset setted ge.0 means that the 
  !     program will encrease rend if needed in such a way that:
  !            [l(i)*(l(i)+1)/rend].lt.[(e-e(i))/irxset ]
  ! rswitch is the switching R between the short range and the 
  ! long range propagator
  ! eps is accuracy for var=ph propagator when the K matrix is calculated
  !     using automatic end of the propagation.
  ! sigch and soglia are accuracy for var-ph propagators: the sigch 
  ! control the channel reduction procedure. 
  ! onlymod is .true. if only modified equations are used in var-ph,
  ! i.e. no automatic end of the propagation is considered. 
  vib_levels=num_vib_levels
  rbegin=0.5_wp
  rswitch=5.0_wp
  rend=10.0_wp
  steps=20
  irxset=0.0_wp
  eps=1e-4
  soglia=1e-6
  sigch=1e-5
  onlymod=.true.
  use_variable_phase=.true.
  care=.false.
 if (myrank==0) then
     read(*, nml = propagator)
  endif

#ifdef PARALLEL
  call mpi_bcast(steps,  1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(irxset, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(rend,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(rswitch,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(rbegin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(eps,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(soglia, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(sigch,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(use_variable_phase,1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(onlymod,           1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(care,              1,MPI_LOGICAL,0, MPI_COMM_WORLD, ierr_mpi)
#endif

  ! +-----------------------------+
  ! |   angular momentum values   |
  ! |   namelist partial_waves    |
  ! +-----------------------------+
  xjtotmin=0.0_wp
  xjtotmax=10.0_wp
  xjtotstep=1.0_wp

  if (myrank==0) then
     read(*, nml = partial_waves)
  endif

! consistency check on the partial waves
  if(xjtotmin>xjtotmax) then
  write(*,*)''
  write(*,*)'Jmin must not be greater than Jmax. Please check it.'
  write(*,*)'Program stop.'
  stop 
  endif

  if(singlet.or.triplet) then
   if((xjtotmin-real(int(xjtotmin)))/=0.0_wp &
      .or.(xjtotmax-real(int(xjtotmax)))/=0.0_wp)then
   write(*,*)''
   write(*,*)'Non-physical choice of Jtotmin or Jtotmax. Please check it.'
   write(*,*)''
   write(*,*)'Program stop.'
   stop
   endif

  endif

  if(doublet) then
   if((xjtotmin-real(int(xjtotmin)))/=0.5_wp &
      .or.(xjtotmax-real(int(xjtotmax)))/=0.5_wp)then
   write(*,*)''
   write(*,*)'Non-physical choice of Jtotmin or Jtotmax. Please check it.'
   write(*,*)''
   write(*,*)'Program stop.'
   stop
   endif
  xjtotmin=xjtotmin+corr
  xjtotmax=xjtotmax+corr
  endif

  if(triplet.and.xjtotmin==0.0_wp) then
   write(*,*)''
   write(*,*)'Jtotmin cannot be 0 in a triplet state of spin.'
   write(*,*)'Please check it.'
   write(*,*)''
   write(*,*)'Program stop.'
   stop
  endif

  if(xjtotstep/=1.0.and.xjtotstep/=2.0)then
   write(*,*)''
   write(*,'(x,a,f3.1)')'xjtotstep=',xjtotstep
   write(*,*)'Strange choice of J step. Please check it.'
   write(*,*)''
   write(*,*)'Program stop.'
   stop
  endif

  jtotmin=nint(xjtotmin)
  jtotmax=nint(xjtotmax)
  jtotstep=nint(xjtotstep)

#ifdef PARALLEL
  call mpi_bcast(jtotmin,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(jtotmax,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(jtotstep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
#endif  

  ! +-----------------------------+
  ! |   extrapolation exponents   |
  ! |   namelist exponents        |
  ! +-----------------------------+
  ! These are the exponents with wich we calculate the correct multipolar
  ! expansion of the potential at long range. You have to supply non-zero 
  ! values for each v_lambda (independently from the n,n' values) you want 
  ! to correctly describe at long range.
  ! The values setted to 0 or not supplied in namelist are setted
  ! directly to "fast" 10 and 12.
  !------------------------------------------------------------------------
  !For example for a neutral_atom-neutral_molecule collision 
  !using 4 v_lambda you can use something similar to:
  !   nlr(1)=6, nlr(2)=6, nlr(3)=6, nlr(4)=0,
  !   mlr(1)=8, mlr(2)=8, mlr(3)=8, mlr(4)=0
  !  (lambda=0, lambda=1, lambda=2, lambda=3)
  !V_0 -> A*R**(-6)+B*R**(-8)
  !V_2 -> A*R**(-6)+B*R**(-8)
  !V_4 -> A*R**(-6)+B*R**(-8)
  !V_6 -> A*R**(-10)+B*R**(-12)
  !   Real potential at long range (with V(4) ->0)
  !   V=C_6(theta)*R**(-6) + C_8(theta)*R**(-8)
  !   where C_6 amd C_8 are assumed to have an anisotropy
  !   rapresentable with the first 3 v_lambda while we let 
  !   the fourth go to 0 very rapidly.  
  !------------------------------------------------------------------------
  !for example for a charged_atom-neutral_molecule collision 
  !using 4 v_lambda you have to use something similar to:
  !   nlr(1)=4, nlr(2)=2, nlr(3)=3, nlr(4)=6,
  !   mlr(1)=6, mlr(2)=6, mlr(3)=4, mlr(4)=8
  !  (lambda=0, lambda=1, lambda=2, lambda=3)
  !V_0 -> A*R**(-4)+B*R**(-6)
  !V_1 -> A*R**(-2)+B*R**(-6)
  !V_2 -> A*R**(-3)+B*R**(-4)
  !V_3 -> A*R**(-6)+B*R**(-8)
  !   Real potential at long range (V(3)->0)
  !   spherical_polarizability + charge  q * a0 * P_0 * R**(-4)
  !   dipole + charge                    q * mu * P_1 * R**(-2)
  !   symmetric_polarizability + charge  q * a2 * P_2 * R**(-4)
  !   quadrupole + charge                q * Q  * P_3 * R**(-3)     
  !   some dispersion in P_0, P_1 and P_3 is added                         
  !   If you have an homonuclear molecule just use only even terms
  !-----------------------------------------------------------------------
  !for example for a neutral_atom-charged_molecule collision 
  !using 4 v_lambda you'll have to use something similar to:
  !   nlr(1)=4, nlr(2)=5, nlr(3)=6, nlr(4)=0,
  !   mlr(1)=6, mlr(2)=6, mlr(3)=8, mlr(4)=0
  !  (lambda=0, lambda=1, lambda=2, lambda=3)
  !V_0 -> A*R**(-4)+B*R**(-6)
  !V_2 -> A*R**(-5)+B*R**(-6)
  !V_4 -> A*R**(-6)+B*R**(-8)
  !V_6 -> A*R**(-10)+B*R**(-12)
  !   Real potential at long range (V(3) and V(4) ->0)
  !   spherical_polarizability + charge  q * a0 * P_0 * R**(-4)
  !   induced dipole + charge            q * mu * P_1 * R**(-5)
  !   some dispersion in P_0, P_1 and P_2 is added
  !   P_4 symmetry is let to go to zero very quickly
  allocate(n_exp(num_lambda_terms), m_exp(num_lambda_terms))
  nlr=0; mlr=0 
  if (myrank==0) then
     read(*, nml= exponents) 
     n_exp(1:num_lambda_terms)=nlr(1:num_lambda_terms)
     m_exp(1:num_lambda_terms)=mlr(1:num_lambda_terms)
     ! for zero or missing value means assume fast decay to zero
     where(n_exp==0) n_exp=10 
     where(m_exp==0) m_exp=12
  endif

#ifdef PARALLEL
  call mpi_bcast(n_exp,  num_lambda_terms, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(m_exp,  num_lambda_terms, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi) 
  !end of the input: syncronize
  call mpi_barrier(MPI_COMM_WORLD, ierr_mpi)
#endif

#ifdef PARALLEL
  write(*,*)''
  write(*,666)'Process', myrank, 'of', nproc, 'is alive: start ...'
  write(*,*)''
666 format(A,1X,I2,1X,A,1X,I2,1X,A)
#endif

  ! debugging unit:
  ! fort.20 rank 0 debugging messages
  ! fort.21 rank 1 debugging messages
  ! fort.22 rank 2 debugging messages
  ! .................................
  unit_mpi=20+myrank
  if (nproc==1) unit_mpi=6  !std output for "dummy mpi"

#ifdef PARALLEL
  ! allocate space for the jtot loops on the various processes
  ! each process will do jtot_per_proc(rank+1) partial waves from 
  ! jtotm(rank+1) to jtotx(rank+1)
  allocate (jtotm(nproc),jtotx(nproc),jtot_per_proc(nproc))  
#endif

  ! allocate quantum numbers and energies vectors for the asymptotic states 
  allocate(e_level(num_N_levels*num_vib_levels))  
  allocate   (n_qn(num_N_levels*num_vib_levels))  
  allocate   (j_qn(num_N_levels*num_vib_levels))  
  allocate   (r_qn(num_N_levels*num_vib_levels))  
  allocate   (s_qn(num_N_levels*num_vib_levels))

  !  allocate several cross sections arrays 
  allocate(symmetry_xs2(num_vib_levels, num_vib_levels, &
   num_rot_levels, num_rot_levels, 3, 3))

  allocate(partial_xs2 (num_vib_levels, num_vib_levels, &
   num_rot_levels, num_rot_levels, 3, 3))
  allocate(total_xs2  (num_vib_levels, num_vib_levels, &
   num_rot_levels, num_rot_levels, 3, 3))

  allocate(sigma_J2(num_vib_levels, num_vib_levels,&
   num_rot_levels, num_rot_levels, 3, 3, 0:jtotmax))
                   
  allocate(sigma_JJ2(0:jtotmax))
  allocate(sigma_el2(num_vib_levels),sigma_in2(num_vib_levels))

  allocate(xxtmp2(num_vib_levels,num_vib_levels))

  ! calculates first 10000 factorial stored in common /fct/
  ! this is needed by xf3j and xf6j.
  call factorial

  ! convert all input into a.u.
  two_mu=2.0_wp*reduced_mass*amu_to_au
  vibr_energy  = (vibr_energy - vibr_energy(1)) * cm_to_au
  rot_constant = rot_constant * cm_to_au
  spin_spin    = spin_spin    * cm_to_au 
  spin_rot     = spin_rot     * cm_to_au 
  rbegin       = rbegin  / au_to_angst
  rend         = rend    / au_to_angst
  rswitch      = rswitch / au_to_angst 
  first_energy = first_energy * cm_to_au
  if(.not.logarithmic_e) delta_e  = delta_e * cm_to_au

  ! here we decide j_step and initial_j in all the cases:

  ! cases: homonuclear + rot_state_symm=0 + doublet
  !        homonuclear + rot_state_symm=0 + triplet
  !        homonuclear + rot_state_symm=1 + triplet
  !        heteronuclear in all cases 
  j_step=1
  initial_j=0

  ! rest of cases:
  if (homonuclear.and.rot_state_symm==0.and.singlet) then
  j_step=2
  initial_j=0
  endif

  if (homonuclear.and.rot_state_symm==1.and.singlet) then
  j_step=2
  initial_j=1
  endif

  if (homonuclear.and.rot_state_symm==0.and.doublet) then
  j_step=1
  initial_j=0
  endif

  ! +---------------------------------+
  ! |    writing out std output       |
  ! | only process 0 is active in i/o |
  ! +---------------------------------+

  rank0_output: if (myrank==0) then
     write(*,1000)'+-----------------------------------------------+ '
     write(*,1000)'|                                               | '
     write(*,1000)'|                    ASPIN                      | '
     write(*,1000)'|                                               | '
     write(*,1000)'| D. Lopez-Duran, E. Bodo, and F. A. Gianturco. | '
     write(*,1000)'|                                               | '
     write(*,1000)'+-----------------------------------------------+ '
     write(*,*)' '
     write(*,*)'Main program started for vibrating diatom-atom &
                  scattering process.'
     write(*,*)'Vibrating diatom in singlet, doublet or triplet &
                  state of spin.'
     write(*,*)'Atom in 1S state. Full close coupling calculations.'
     write(*,*)' '
     write(*,1000)'+-----------------------------------+'
     write(*,1000)'|   namelist: molecular_general     |'
     write(*,1000)'+-----------------------------------+'
     write(*,*)' '
     write(*,1001)' Vib. levels=',num_vib_levels, &
                  'Max j value=',xmax_j_value, &
                  'Initial vib. state=',initial_vib_level-1
     write(*,*)' '
     write(*,998)' Initial j state=', xinit_j_state, &
                 'Initial N state=', init_N_state
     write(*,*)''
     write(*,994)'Reduced mass of the system (amu) =',reduced_mass
     write(*,*)''
     write(*,1034)'Spin state of the molecule = ', stsp
     write(*,*)''
     write(*,*)'(relaxed=',relaxed,')'
     write(*,*)''
     write(*,1000)'+-----------------------------------+'
     write(*,1000)'|        namelist: energy           |'
     write(*,1000)'+-----------------------------------+'
     write(*,*)''
     write(*,995)'First energy (cm-1) =',first_energy / cm_to_au
     write(*,*)''
     write(*,995)'Delta_energy (cm-1) =',delta_e / cm_to_au
     write(*,*)''
     write(*,992)'Number of energies=',num_of_energies
     write(*,*)''
     write(*,*)'(logarithmic_e=',logarithmic_e,')'
     write(*,*)''
     write(*,1000)'+-----------------------------------+'
     write(*,1000)'|   namelist: potential_general     |'
     write(*,1000)'+-----------------------------------+'
     write(*,*)''
     write(*,992)'Number of lambda terms =',num_lambda_terms
     write(*,*)''
     write(*,992)'num_pot_points=',num_pot_points
     write(*,*)''
     write(*,992)'Points to discard=',points_to_discard
     write(*,*)''
     write(*,992)'lambda_to_zero=',lambda_to_zero
     write(*,*)''
     write(*,1000)'+-----------------------------------+'
     write(*,1000)'|   namelist: molecular_params      |'
     write(*,1000)'+-----------------------------------+'
     write(*,*)''
     write(*,*)'Vibrational energies:'
     write(*,1006)'v','e(v)/cm-1'
     do i=1,num_vib_levels
        write(*,1007)i-1,vibr_energy(i)/cm_to_au
     enddo
     write(*,*)''
     write(*,*)'Rotational constants:'
     write(*,1002)'v','b(v)/cm-1'
     do i=1,num_vib_levels
        write(*,1003)i-1,rot_constant(i)/cm_to_au
     enddo
     write(*,*)''
     write(*,2002)'Spin-spin interaction:','Lambda/cm-1  '
     write(*,2003) spin_spin/cm_to_au
     write(*,*)' '
     write(*,2002)'Spin-rotation interaction:','Gamma/cm-1  '
     write(*,2003) spin_rot/cm_to_au
     write(*,*)''
     write(*,*)'(centrifugal distortion is neglected)'
     write(*,*)''
     write(*,1004)'For level:','v','max j value is'
     do i=1,num_vib_levels
        write(*,1005)i-1,xjmax(i)
     enddo
     write(*,*)''
     write(*,*)'(homonuclear=',homonuclear,')'
     write(*,*)''
     write(*,993)'(rot_state_symm=',rot_state_symm,')'
     write(*,*)''
     write(*,1000)'+-----------------------------------+'
     write(*,1000)'|   namelist: program_behaviour     |'
     write(*,1000)'+-----------------------------------+'
     write(*,*)''
     write(*,992)'Level of debugging=',debug
     write(*,*)''
     write(*,992)'Unit for the potential=',potential_unit
     write(*,*)''
     write(*,1000)'+-----------------------------------+'
     write(*,1000)'|       namelist: propagator        |'
     write(*,1000)'+-----------------------------------+'
     write(*,*)''
     write(*,991)'rbegin (angstrom) = ', rbegin * au_to_angst
     write(*,*)''
     write(*,991)'rswitch (angstrom) = ', rswitch * au_to_angst
     write(*,*)''
     write(*,991)'rend (angstrom) = ', rend * au_to_angst
     write(*,*)''
     write(*,992)'Number of steps = ', steps
     write(*,*)''
     write(*,991)'irxset parameter = ', irxset
     write(*,*)''
     write(*,991)'epsilon parameter =',eps
     write(*,*)''
     write(*,991)'sigma parameter =',sigch
     write(*,*)''
     write(*,*)'(use_variable_phase =',use_variable_phase,')'
     write(*,*)''
     write(*,*)'(onlymod =',onlymod,')'
     write(*,*)''
     write(*,*)'(care =',care,')'
     write(*,*)''
     write(*,1000)'+-----------------------------------+'
     write(*,1000)'|     namelist: partial_waves       |'
     write(*,1000)'+-----------------------------------+'
     write(*,*)''
     write(*,996)'Minimum total J =',xjtotmin-corr
     write(*,*)''
     write(*,996)'Maximum total J =',xjtotmax-corr
     write(*,*)''
     write(*,996)'Step in J =',xjtotstep
     write(*,*)''
     write(*,1000)'+-----------------------------------+'
     write(*,1000)'|       namelist: exponents         |'
     write(*,1000)'+-----------------------------------+'
     write(*,*)''
     write(*,*)'Exponents:'
     write(*,1046)'i','nlr(i)','mlr(i)'
     do i=1,num_lambda_terms
        write(*,1047)i,nlr(i),mlr(i)
     enddo

     ! consistency check on the propagator variables 
     ! we assume use_variable_phase in strange case
     if (rbegin>=rswitch) then
        write(*,*)'CRITICAL: rbegin > rswitch'
        write(*,*)'Program stopping'

#ifdef PARALLEL
        write(*,*)'stopping mpi subsystem'
        call mpi_finalize(ierr_mpi)
#endif
        stop
     endif

     ! rend will be taken into consideration only
     ! if use_variable_phase is true 
     if (use_variable_phase) then
        if (rbegin>=rend) then
           write(*,*)'CRITICAL: rbegin > rend'
           write(*,*)'Program stopping'

#ifdef PARALLEL
           write(*,*)'stopping mpi subsystem'
           call mpi_finalize(ierr_mpi)
#endif
           stop
        endif
        if (rbegin>=rswitch) then
           write(*,*)'CRITICAL: rbegin > rswitch'
           write(*,*)'Program stopping'

#ifdef PARALLEL
           write(*,*)'stopping mpi subsystem'
           call mpi_finalize(ierr_mpi)
#endif
           stop
        endif
        if (rswitch>=rend) then
           write(*,*)'WARNING: rswitch > rend'
           write(*,*)'Propagation from rswitch to rend &
                      will not be performed'
           write(*,*)'Only pure log-derivative propagator will be setted'
           use_variable_phase=.false.
        endif
     endif

     write(*,*)' '
     write(*,*)' '
     write(*,1000)'+-----------------------------------------------+'
     write(*,1000)'|propagator and potential parameter(in angtroms)|'
     write(*,1000)'+-----------------------------------------------+'
     write(*,*)' '
     write(*,1008)'number of vibrational couplings','=', & 
          num_vib_levels*(num_vib_levels+1)/2
     write(*,1008)'number of points','=',num_pot_points
     write(*,1008)'number of terms in the lambda expansion', &
                  '=',num_lambda_terms
     write(*,*)' '
     write(*,1009)'initial point for integration','=',rbegin*au_to_angst
     write(*,1009)'matching point requested','=',rend*au_to_angst

     if (use_variable_phase) then
        write(*,fmt='(/,t10,a,/)')'Hibrid Manolopulos Var-phase &
                                   propagator has been chosen'
        write(*,1009)'switching point between PM and MA','=',&
             rswitch*au_to_angst
     else
        write(*,fmt='(/,t10,a,/)')'Pure Manolopulos propagator &
                                   has been chosen'
     endif

     write(*,1008)'step parameter for log-der propagator','=',steps

     if (points_to_discard/=0) then
        write(*,999)'warning: points_to_discard !=0'
        write(*,999)'you requested to use inner points in r to get'
        write(*,999)'    the long range fitting'
     endif

897 format(x,a,i3,2x,a,f5.1)
898 format(x,a,f4.1,5x,a,f4.1)
899 format(x,a,f5.1)
990 format(x,a,f4.1,5x,a,i3)
991 format(x,a,f11.5)
992 format(x,a,i6)
993 format(x,a,i2,a)
994 format(x,a,f7.5,a)
995 format(a22,es11.4)
996 format(x,a,f4.1)
997 format(a9,2x,f4.1)
998 format(a17,f4.1,5x,a16,i3)
999  format(t1,a)
1000 format(t3,a)
1001 format(a12,i3,4x,a12,f4.1,4x,a19,1x,i1)
1002 format(t25,a2,t30,a9)
2002 format(x,a,2x,t30,a13)
2003 format(t30,es11.4)
1003 format(t25,i2,t30,f9.4)
1004 format(x,a,2x,t25,a2,t30,a15)
1005 format(t25,i2,t30,f9.1)
1006 format(t25,a2,t30,a9)
1007 format(t25,i2,t30,f9.2)
1008 format(a,t43,a1,t45,i6)
1009 format(a,t43,a1,t45,f15.4)
1030 format(a,1x,f3.1)
1033 format(a,i2)
1034 format(a,a)
1046 format(t21,a1,t26,a6,t34,a6)
1047 format(t20,i2,t27,i2,t35,i2)

  endif rank0_output
  
  ! +-------------------------------------------+
  ! |   mpi: control return to all processes    |
  ! |    ordering asymptotic channels:          |
  ! +-------------------------------------------+
  count=0                    
  max_j_for_n = max_j_for_n + 1   
  initial_j   = initial_j   + 1       
  ! max_j_for_n <0 means to use max_j_value in each vibrational manifold
  ! we stop if we have an incoherent max_j_for_n
  do i=1,num_vib_levels
      if (max_j_for_n(i)>num_rot_levels) then
        write(*,*)'error in main: max_j_for_n must be less than max_j_value'

#ifdef PARALLEL
        write(*,*)'stopping mpi subsystem'
        call mpi_finalize(ierr_mpi)
#endif

        stop
     endif
     if (max_j_for_n(i)<=0) then
        write(*,1010)'for vib. level',i-1,'max_j_for_n is .lt. zero'
        write(*,1010)'we set',max_j_value-1,'as the maximum allowed j'
        max_j_for_n(i)=num_rot_levels
     endif
  enddo

  ! start cycling - counting channels
  if (myrank==0)then
     write(*,*)' '
     write(*,1000)'+------------------------------+'
     write(*,1000)'| main:  **ordering channels** |'
     write(*,1000)'+------------------------------+'
     write(*,*)' '
  endif

  if (relaxed) then
     if (myrank==0)then
        write(*,999)'main: flag relaxed=.true.'
        write(*,999)'main: warning using the ro-vibrational &
                     spectrum as it is'
        write(*,999)'      nearly degenerate channels could be troublesome'
     endif
  endif

  vib_loop: do n=1,num_vib_levels
     if (myrank==0)then
        write(*,*)' '
        write(*,1011)'vib level no','|','min_j','|','max_j','|','j_step'
        write(*,*)'+------------------------------------------+'
        write(*,1012)n-1,'|',initial_j-1,'|', & 
                          real(max_j_for_n(n)-1, kind=wp)-corr,'|',j_step
        write(*,*)'+------------------------------------------+'
        write(*,1013)'index','|','v','|','j','|','N','|','energy/cm-1'
     endif

     uno=1.0_wp
     dos=2.0_wp
     tres=3.0_wp

     j_loop: do j=initial_j, max_j_for_n(n), j_step

! HOMONUCLEAR CASE 

      homo: if (homonuclear) then  

! HOMONUCLEAR CASE ANR ROT_STATE_SYMM = 0

             if (rot_state_symm==0) then ! homonuclear + rot state = 0

             three_cases_a: if(singlet)then
             !singlet case
             jreal=j-1

         e=vibr_energy(n)+rot_constant(n)*real(jreal*(jreal+1),kind=wp)
         ! check that we are not going on the next vib level with
         ! rotational channels if relaxed=.false. (default)
         ! we use all the j in the last vib. level anyway
         if (.not.relaxed) then
          if(n<num_vib_levels)  e_to_compare=vibr_energy(n+1)
          if(n==num_vib_levels) e_to_compare=1.0e+100_wp
          if (e<e_to_compare) then
           count=count+1
           e_level(count)=e      ! filling energy vectors
           n_qn(count)=n-1       ! filling vib. quantum numbers vectors
           j_qn(count)=jreal     ! filling rot. quantum numbers vectors
           r_qn(count)=jreal     ! filling N quantum numbers vectors
           s_qn(count)=2
            if (myrank==0)then
             write(*,1014)count,'|',n_qn(count),'|', &
             real(j_qn(count),kind=wp),'|',r_qn(count), &
              '|',e_level(count)/cm_to_au
            endif
          else
           if (myrank==0)then
            write(*,1015)'*','|',jreal,'|',n-1,'|', &
                         jreal,'|',e/cm_to_au 
           endif
          endif
         else  !taking the spectrum as it is
          count=count+1
          e_level(count)=e      ! filling energy vectors
          n_qn(count)=n-1       ! filling vib. quantum numbers vectors
          j_qn(count)=jreal     ! filling rot. quantum numbers vectors
          r_qn(count)=jreal     ! filling N quantum numbers vectors
          s_qn(count)=2
          if (myrank==0)then
           write(*,1014)count,'|',n_qn(count),'|', &
           real(j_qn(count),kind=wp),'|',r_qn(count), &
           '|',e_level(count)/cm_to_au
          endif
         endif

! we are in homonuclear + rot state = 0

             elseif(doublet)then
             !doublet case
             jreal=j
              do kr=j,j+1 ! doublet 
               if (parity(kr)) cycle 
                rreal=kr-1
                e=vibr_energy(n)
                 if (rreal<jreal)then
                  ! F1  N=j-1/2  
                  e=e+rot_constant(n)*rreal*(rreal+uno)
                  e=e+.5d0*spin_rot*rreal
                 elseif (rreal==jreal)then
                  ! F2  N=j+1/2 
                  e=e+rot_constant(n)*rreal*(rreal+uno)
                  e=e-.5d0*spin_rot*(rreal+uno)
                 endif
         ! check that we are not going on the next vib level with
         ! rotational channels if relaxed=.false. (default)
         ! we use all the j in the last vib. level anyway
                 if (.not.relaxed) then
                  if(n<num_vib_levels)  e_to_compare=vibr_energy(n+1)
                  if(n==num_vib_levels) e_to_compare=1.0e+100_wp
                  if (e<e_to_compare) then
                   count=count+1
                   e_level(count)=e  ! filling energy vectors
                   n_qn(count)=n-1   ! filling vib. q. numbers vectors
                   j_qn(count)=jreal ! filling rot. q. numbers vectors
                   r_qn(count)=rreal ! filling N q. numbers vectors
                   if(jreal==rreal)then
                   s_qn(count)=2
                   elseif(jreal>rreal)then
                   s_qn(count)=3
                   endif
                    if (myrank==0)then
                     write(*,1014)count,'|',n_qn(count),'|', &
                     real(j_qn(count),kind=wp)-corr,'|', &
                     r_qn(count), '|',e_level(count)/cm_to_au
                    endif
                   else
                    if (myrank==0)then
                    write(*,1015)'*','|',jreal,'|',n-1, &
                                 '|',rreal,'|',e/cm_to_au
                    endif
                   endif
                  else  !taking the spectrum as it is
                   count=count+1
                   e_level(count)=e  ! filling energy vectors
                   n_qn(count)=n-1   ! filling vib. q. numbers vectors
                   j_qn(count)=jreal ! filling rot. q. numbers vectors
                   r_qn(count)=rreal ! filling N q. numbers vectors
                   if(jreal==rreal)then
                   s_qn(count)=2
                   elseif(jreal>rreal)then
                   s_qn(count)=3
                   endif
                   if (myrank==0)then
                    write(*,1014)count,'|',n_qn(count),'|', &
                    real(j_qn(count),kind=wp)-corr,'|',r_qn(count), &
                    '|',e_level(count)/cm_to_au
                   endif
                  endif
              enddo ! doublet

! we are in homonuclear + rot state = 0

             elseif(triplet)then
             !triplet case
             jreal=j-1
              do kr=j,j+2 ! triplet 
               if (.not.parity(kr)) cycle
               if (j==1.and.kr<3) goto 13 
               rreal=kr-2
               e=vibr_energy(n)
                if (rreal<jreal)then
                 ! F1  N=j-1
                 e=e+rot_constant(n)*rreal*(rreal+uno)
                 e=e+rot_constant(n)*(dos*rreal+tres)
                 e=e-spin_spin
                 e=e-sqrt(rot_constant(n)**2*(dos*rreal+tres)**2 &
                          +spin_spin**2-dos*spin_spin*rot_constant(n))
                 e=e+spin_rot*(rreal+uno)
                elseif (rreal==jreal)then
                 ! F2  N=j
                 e=e+rot_constant(n)*rreal*(rreal+uno)
                elseif (rreal>jreal)then
                 ! F3  N=j+1
                 e=e+rot_constant(n)*rreal*(rreal+uno)
                 e=e-rot_constant(n)*(dos*rreal-uno)
                 e=e-spin_spin
                 e=e+sqrt(rot_constant(n)**2*(dos*rreal-uno)**2 &
                          +spin_spin**2-dos*spin_spin*rot_constant(n))
                 e=e-spin_rot*rreal
                endif
         ! check that we are not going on the next vib level with
         ! rotational channels if relaxed=.false. (default)
         ! we use all the j in the last vib. level anyway
                if (.not.relaxed) then
                if(n<num_vib_levels)  e_to_compare=vibr_energy(n+1)
                if(n==num_vib_levels) e_to_compare=1.0e+100_wp
                if (e<e_to_compare) then
                 count=count+1
                 e_level(count)=e    ! filling energy vectors
                 n_qn(count)=n-1     ! filling vib. q. numbers vectors
                 j_qn(count)=jreal   ! filling rot. q. numbers vectors
                 r_qn(count)=rreal   ! filling N q. numbers vectors
                 if(jreal<rreal)then
                 s_qn(count)=1
                 elseif(jreal==rreal)then
                 s_qn(count)=2
                 elseif(jreal>rreal)then
                 s_qn(count)=3
                 endif
                  if (myrank==0)then
                   write(*,1014)count,'|',n_qn(count),'|', &
                   real(j_qn(count),kind=wp),'|',r_qn(count), &
                   '|',e_level(count)/cm_to_au
                  endif
                else
                if (myrank==0)then
                   write(*,1015)'*','|',jreal,'|',n-1,'|', &
                   rreal,'|',e/cm_to_au
                endif
               endif
              else  !taking the spectrum as it is
               count=count+1
                e_level(count)=e      ! filling energy vectors
                n_qn(count)=n-1       ! filling vib. q. numbers vectors
                j_qn(count)=jreal     ! filling rot. q. numbers vectors
                r_qn(count)=rreal     ! filling N q. numbers vectors
                 if(jreal<rreal)then
                 s_qn(count)=1
                 elseif(jreal==rreal)then
                 s_qn(count)=2
                 elseif(jreal>rreal)then
                 s_qn(count)=3
                 endif
                 if (myrank==0)then
                  write(*,1014)count,'|',n_qn(count),'|', &
                  real(j_qn(count),kind=wp),'|',r_qn(count), &
                  '|',e_level(count)/cm_to_au
                 endif
              endif

13 continue

              enddo !triplet

         endif three_cases_a

              ! end of homonuclear + rot state = 0

! HOMONUCLEAR CASE AND ROT_STATE_SYMM = 1
 
             else if (rot_state_symm==1) then ! homonuclear + rot state = 1

         three_cases_b: if(singlet)then
         !singlet case
         jreal=j-1

         e=vibr_energy(n)+rot_constant(n)*real(jreal*(jreal+1),kind=wp)
         ! check that we are not going on the next vib level with
         ! rotational channels if relaxed=.false. (default)
         ! we use all the j in the last vib. level anyway
         if (.not.relaxed) then
          if(n<num_vib_levels)  e_to_compare=vibr_energy(n+1)
          if(n==num_vib_levels) e_to_compare=1.0e+100_wp
          if (e<e_to_compare) then
           count=count+1
           e_level(count)=e      ! filling energy vectors
           n_qn(count)=n-1       ! filling vib. q. numbers vectors
           j_qn(count)=jreal     ! filling rot. q. numbers vectors
           r_qn(count)=jreal     ! filling N q. numbers vectors
           s_qn(count)=2
           if (myrank==0)then
            write(*,1014)count,'|',n_qn(count),'|', &
            real(j_qn(count),kind=wp),'|',r_qn(count), &
            '|',e_level(count)/cm_to_au
           endif
          else
           if (myrank==0)then
            write(*,1015)'*','|',jreal,'|',n-1,'|', &
                         jreal,'|',e/cm_to_au
           endif
          endif
         else  !taking the spectrum as it is
          count=count+1
          e_level(count)=e      ! filling energy vectors
          n_qn(count)=n-1       ! filling vib. q. numbers vectors
          j_qn(count)=jreal     ! filling rot. q. numbers vectors
          r_qn(count)=jreal     ! filling N q. numbers vectors
          s_qn(count)=2
          if (myrank==0)then
           write(*,1014)count,'|',n_qn(count),'|', &
           real(j_qn(count),kind=wp),'|',r_qn(count), &
           '|',e_level(count)/cm_to_au
           endif
         endif

! we are in homonuclear + rot state = 1

             elseif(doublet)then
             !doublet case
              jreal=j
              do kr=j,j+1 ! doublet
               if (.not.parity(kr)) cycle
               rreal=kr-1
               e=vibr_energy(n)
                if (rreal<jreal)then
                 ! F1  N=j-1/2  
                 e=e+rot_constant(n)*rreal*(rreal+uno)
                 e=e+.5d0*spin_rot*rreal
                elseif (rreal==jreal)then
                 ! F2  N=j+1/2  
                 e=e+rot_constant(n)*rreal*(rreal+uno)
                 e=e-.5d0*spin_rot*(rreal+uno)
                endif
         ! check that we are not going on the next vib level with
         ! rotational channels if relaxed=.false. (default)
         ! we use all the j in the last vib. level anyway
                if (.not.relaxed) then
                if(n<num_vib_levels)  e_to_compare=vibr_energy(n+1)
                if(n==num_vib_levels) e_to_compare=1.0e+100_wp
                 if (e<e_to_compare) then
                  count=count+1
                  e_level(count)=e   ! filling energy vectors
                  n_qn(count)=n-1    ! filling vib. q. numbers vectors
                  j_qn(count)=jreal  ! filling rot. q. numbers vectors
                  r_qn(count)=rreal  ! filling N q. numbers vectors
                  if(jreal==rreal)then
                  s_qn(count)=2
                  elseif(jreal>rreal)then
                  s_qn(count)=3
                  endif
                   if (myrank==0)then
                    write(*,1014)count,'|',n_qn(count),'|', &
                    real(j_qn(count),kind=wp)-corr,'|',r_qn(count), &
                    '|',e_level(count)/cm_to_au
                   endif
                 else
                  if (myrank==0)then
                   write(*,1015)'*','|',jreal,'|',n-1,'|', &
                                rreal,'|',e/cm_to_au
                  endif
                 endif
                else  !taking the spectrum as it is
                 count=count+1
                 e_level(count)=e    ! filling energy vectors
                 n_qn(count)=n-1     ! filling vib. q. numbers vectors
                 j_qn(count)=jreal   ! filling rot. q. numbers vectors
                 r_qn(count)=rreal   ! filling N q. numbers vectors
                 if(jreal==rreal)then
                 s_qn(count)=2
                 elseif(jreal>rreal)then
                 s_qn(count)=3
                 endif
                 if (myrank==0)then
                  write(*,1014)count,'|',n_qn(count),'|', &
                  real(j_qn(count),kind=wp)-corr,'|',r_qn(count), &
                  '|',e_level(count)/cm_to_au
                 endif
                endif
             enddo ! doublet

! we are in homonuclear + rot state = 1

             elseif(triplet)then
             !triplet case
             jreal=j-1
             do kr=j,j+2 ! triplet
              if (parity(kr)) cycle
              if (j==1.and.kr<3) goto 14
              rreal=kr-2
              e=vibr_energy(n)
               if (rreal<jreal)then
                ! F1  N=j-1
                e=e+rot_constant(n)*rreal*(rreal+uno)
                e=e+rot_constant(n)*(dos*rreal+tres)
                e=e-spin_spin
                e=e-sqrt(rot_constant(n)**2*(dos*rreal+tres)**2 &
                         +spin_spin**2-dos*spin_spin*rot_constant(n))
                e=e+spin_rot*(rreal+uno)
               elseif (rreal==jreal)then
                ! F2  N=j
                e=e+rot_constant(n)*rreal*(rreal+uno)
               elseif (rreal>jreal)then
                ! F3  N=j+1
                e=e+rot_constant(n)*rreal*(rreal+uno)
                e=e-rot_constant(n)*(dos*rreal-uno)
                e=e-spin_spin
                e=e+sqrt(rot_constant(n)**2*(dos*rreal-uno)**2 &
                         +spin_spin**2-dos*spin_spin*rot_constant(n))
                e=e-spin_rot*rreal
               endif
          ! check that we are not going on the next vib level with
          ! rotational channels if relaxed=.false. (default)
          ! we use all the j in the last vib. level anyway
               if (.not.relaxed) then
               if(n<num_vib_levels)  e_to_compare=vibr_energy(n+1)
               if(n==num_vib_levels) e_to_compare=1.0e+100_wp
               if (e<e_to_compare) then
                count=count+1
                e_level(count)=e      ! filling energy vectors
                n_qn(count)=n-1       ! filling vib. q. numbers vectors
                j_qn(count)=jreal     ! filling rot. q. numbers vectors
                r_qn(count)=rreal     ! filling N q. numbers vectors
                 if(jreal<rreal)then
                 s_qn(count)=1
                 elseif(jreal==rreal)then
                 s_qn(count)=2
                 elseif(jreal>rreal)then
                 s_qn(count)=3
                 endif
                if (myrank==0)then
                 write(*,1014)count,'|',n_qn(count),'|', &
                 real(j_qn(count),kind=wp),'|',r_qn(count), &
                 '|',e_level(count)/cm_to_au
                endif
               else
                if (myrank==0)then
                 write(*,1015)'*','|',jreal,'|',n-1,'|', &
                              rreal,'|',e/cm_to_au
                endif
              endif
           else  !taking the spectrum as it is
            count=count+1
            e_level(count)=e      ! filling energy vectors
            n_qn(count)=n-1       ! filling vib. q. numbers vectors
            j_qn(count)=jreal     ! filling rot. q. numbers vectors
            r_qn(count)=rreal     ! filling N q. numbers vectors
             if(jreal<rreal)then
             s_qn(count)=1
             elseif(jreal==rreal)then
             s_qn(count)=2
             elseif(jreal>rreal)then
             s_qn(count)=3
             endif
            if (myrank==0)then
             write(*,1014)count,'|',n_qn(count),'|', &
             real(j_qn(count),kind=wp),'|',r_qn(count), &
             '|',e_level(count)/cm_to_au
            endif

          endif

14 continue

        enddo ! triplet

         endif three_cases_b

           endif ! end of homonuclear + rot state = 1

! HETERONUCLEAR CASE

             else  

             three_cases_c: if(singlet)then
             !singlet case
             jreal=j-1
        e=vibr_energy(n)+rot_constant(n)*real(jreal*(jreal+1),kind=wp)
        ! check that we are not going on the next vib level with
        ! rotational channels if relaxed=.false. (default)
        ! we use all the j in the last vib. level anyway
        if (.not.relaxed) then
         if(n<num_vib_levels)  e_to_compare=vibr_energy(n+1)
         if(n==num_vib_levels) e_to_compare=1.0e+100_wp
         if (e<e_to_compare) then
          count=count+1
          e_level(count)=e      ! filling energy vectors
          n_qn(count)=n-1       ! filling vib. q. numbers vectors
          j_qn(count)=jreal     ! filling rot. q. numbers vectors
          r_qn(count)=jreal     ! filling N q. numbers vectors
          s_qn(count)=2
           if (myrank==0)then
            write(*,1014)count,'|',n_qn(count),'|', &
            real(j_qn(count),kind=wp),'|',r_qn(count), &
            '|',e_level(count)/cm_to_au
            endif
           else
            if (myrank==0)then
             write(*,1015)'*','|',jreal,'|',n-1,'|', &
                          jreal,'|',e/cm_to_au
             endif
           endif
         else  !taking the spectrum as it is
          count=count+1
          e_level(count)=e      ! filling energy vectors
          n_qn(count)=n-1       ! filling vib. q. numbers vectors
          j_qn(count)=jreal     ! filling rot. q. numbers vectors
          r_qn(count)=jreal     ! filling N q. numbers vectors
          s_qn(count)=2
           if (myrank==0)then
            write(*,1014)count,'|',n_qn(count),'|', &
            real(j_qn(count),kind=wp),'|',r_qn(count), &
            '|',e_level(count)/cm_to_au
           endif
        endif

! we are in the heteronuclear case 

             elseif(doublet)then
             !doublet case
             jreal=j
             do kr=j,j+1 ! doublet
              rreal=kr-1
              e=vibr_energy(n)
              if (rreal<jreal)then
               ! F1  N=j-1/2  
               e=e+rot_constant(n)*rreal*(rreal+uno)
               e=e+.5d0*spin_rot*rreal
              elseif (rreal==jreal)then
               ! F2  N=j+1/2  
               e=e+rot_constant(n)*rreal*(rreal+uno)
               e=e-.5d0*spin_rot*(rreal+uno)
              endif
          ! check that we are not going on the next vib level with
          ! rotational channels if relaxed=.false. (default)
          ! we use all the j in the last vib. level anyway
              if (.not.relaxed) then
               if(n<num_vib_levels)  e_to_compare=vibr_energy(n+1)
               if(n==num_vib_levels) e_to_compare=1.0e+100_wp
               if (e<e_to_compare) then
                count=count+1
                e_level(count)=e    ! filling energy vectors
                n_qn(count)=n-1     ! filling vib. q. numbers vectors
                j_qn(count)=jreal   ! filling rot. q. numbers vectors
                r_qn(count)=rreal   ! filling N q. numbers vectors
                 if(jreal==rreal)then
                 s_qn(count)=2
                 elseif(jreal>rreal)then
                 s_qn(count)=3
                 endif 
                if (myrank==0)then
                 write(*,1014)count,'|',n_qn(count),'|', &
                 real(j_qn(count),kind=wp)-corr,'|',r_qn(count), &
                 '|',e_level(count)/cm_to_au
                 endif
                else
                 if (myrank==0)then
                 write(*,1015)'*','|',jreal,'|',n-1,'|', &
                               rreal,'|',e/cm_to_au
                 endif
                endif
              else  !taking the spectrum as it is
               count=count+1
               e_level(count)=e    ! filling energy vectors
               n_qn(count)=n-1     ! filling vib. q. numbers vectors
               j_qn(count)=jreal   ! filling rot. q. numbers vectors
               r_qn(count)=rreal   ! filling N q. numbers vectors
                 if(jreal==rreal)then
                 s_qn(count)=2
                 elseif(jreal>rreal)then
                 s_qn(count)=3
                 endif
               if (myrank==0)then
                write(*,1014)count,'|',n_qn(count),'|', &
                real(j_qn(count),kind=wp)-corr,'|',r_qn(count), &
                '|',e_level(count)/cm_to_au
                endif
              endif  
           enddo ! doublet

! we are in the heteronuclear case

             elseif(triplet)then
             !triplet      
             jreal=j-1     
             do kr=j,j+2 ! triplet
              if (j==1.and.kr<3) cycle
               rreal=kr-2
               e=vibr_energy(n)
               if (rreal<jreal)then
                ! F1  N=j-1
                e=e+rot_constant(n)*rreal*(rreal+uno)
                e=e+rot_constant(n)*(dos*rreal+tres)
                e=e-spin_spin
                e=e-sqrt(rot_constant(n)**2*(dos*rreal+tres)**2 &
                         +spin_spin**2-dos*spin_spin*rot_constant(n))
                e=e+spin_rot*(rreal+uno)
               elseif (rreal==jreal)then
                ! F2  N=j
                e=e+rot_constant(n)*rreal*(rreal+uno)
               elseif (rreal>jreal)then
                ! F3  N=j+1
                e=e+rot_constant(n)*rreal*(rreal+uno)
                e=e-rot_constant(n)*(dos*rreal-uno)
                e=e-spin_spin
                e=e+sqrt(rot_constant(n)**2*(dos*rreal-uno)**2 &
                         +spin_spin**2-dos*spin_spin*rot_constant(n))
                e=e-spin_rot*rreal
               endif
         ! check that we are not going on the next vib level with
         ! rotational channels if relaxed=.false. (default)
         ! we use all the j in the last vib. level anyway
             if (.not.relaxed) then
             if(n<num_vib_levels)  e_to_compare=vibr_energy(n+1)
             if(n==num_vib_levels) e_to_compare=1.0e+100_wp
             if (e<e_to_compare) then
              count=count+1
              e_level(count)=e      ! filling energy vectors
              n_qn(count)=n-1       ! filling vib. q. numbers vectors
              j_qn(count)=jreal     ! filling rot. q. numbers vectors
              r_qn(count)=rreal     ! filling N q. numbers vectors
                 if(jreal<rreal)then
                 s_qn(count)=1
                 elseif(jreal==rreal)then
                 s_qn(count)=2
                 elseif(jreal>rreal)then
                 s_qn(count)=3
                 endif
               if (myrank==0)then
                write(*,1014)count,'|',n_qn(count),'|', &
                real(j_qn(count),kind=wp),'|',r_qn(count), &
                '|',e_level(count)/cm_to_au
                endif
             else
              if (myrank==0)then
               write(*,1015)'*','|',jreal,'|',n-1,'|', &
                            rreal,'|',e/cm_to_au
              endif
             endif
          else  !taking the spectrum as it is
             count=count+1
             e_level(count)=e      ! filling energy vectors
             n_qn(count)=n-1       ! filling vib. q. numbers vectors
             j_qn(count)=jreal     ! filling rot. q. numbers vectors
             r_qn(count)=rreal     ! filling N q. numbers vectors
                if(jreal<rreal)then
                s_qn(count)=1
                elseif(jreal==rreal)then
                s_qn(count)=2
                elseif(jreal>rreal)then
                s_qn(count)=3
                endif
             if (myrank==0)then
              write(*,1014)count,'|',n_qn(count),'|', &
              real(j_qn(count),kind=wp),'|',r_qn(count), &
              '|',e_level(count)/cm_to_au
             endif
          endif
        enddo ! triplet

! end of the heteronuclear case

         endif three_cases_c
        endif homo
     enddo j_loop 
  enddo vib_loop
 
 ! total number of channels assigned to num_ch
  num_ch=count

#ifdef PARALLEL
  ! this is needed to collect data out of the parallel cycle
  if (myrank==0) then
     allocate(real_total2(num_vib_levels, num_vib_levels, &
              num_rot_levels, num_rot_levels, &
              3, 3))
  endif
#endif

  if (myrank==0)then
     write(*,*)' '
     write(*,1016)'main: total number of asymptotic channels:', num_ch
  endif
 
  !--------------------------Bubblesort-------------------------------
  ! we have to sort channels in case relaxed=.true.
 
  if (relaxed) then
     if (myrank==0)then
        write(*,999)'relaxed=.true.' 
        write(*,999)'sorting the energy spectrum'
     endif
     i=0
     do 
        i=i+1
        if(e_level(i).gt.e_level(i+1)) then
           e_new=e_level(i)
           j_new=j_qn(i)
           n_new=n_qn(i)
           r_new=r_qn(i)
           s_new=s_qn(i)
           e_level(i)=e_level(i+1)
           j_qn(i)=j_qn(i+1)
           n_qn(i)=n_qn(i+1)
           r_qn(i)=r_qn(i+1)
           s_qn(i)=s_qn(i+1)
           e_level(i+1)=e_new
           j_qn(i+1)=j_new
           n_qn(i+1)=n_new
           r_qn(i+1)=r_new
           s_qn(i+1)=s_new
           i=0
           cycle
        elseif (i==num_ch-1) then
           exit
        endif
     enddo

     if (myrank==0)then
        write(*,*)' '
        write(*,*)'+------------------------------------------+'
        write(*,1013)'index','|','v','|','j','|','N','|','energy/cm-1'
     endif
     do count=1,num_ch
        if (myrank==0)then
           write(*,1014)count,'|',n_qn(count),'|', &
           real(j_qn(count),kind=wp)-corr,'|',r_qn(count), &
           '|',e_level(count)/cm_to_au
         endif
     enddo
  endif
  
  if (myrank==0)then
     write(*,*)' '
     write(*,999)'main:  ordering channels finished'
     write(*,999)'       note: * labelled channels will not be used'
     write(*,999)'       their energy is higher than vibr_energy(n+1)'
  endif
  !-------------------------------------------------------------------

1010 format(a,t15,i4,t20,a)
1011 format(a,t15,a,2x,a,t25,a,2x,a,t35,a,2x,a)
1012 format(i4,t15,a,i4,t25,a,f6.1,t35,a,i4)
1013 format(a,t15,a,3x,a,t25,a,3x,a,t35,a,3x,a,t45,a,1x,a)
1014 format(i4,t15,a,i4,t25,a,2x,f4.1,t35,a,i4,t45,a,f8.2)
1015 format(a,t15,a,i4,t25,a,i4,t35,a,i4,t45,a,f8.2)
1016 format(a,2x,i4)
1040 format(a,2x,f4.1)
  ! +------------------------------------------------------------+
  ! |   potential treatment and extrapolation, if needed         |  
  ! |                                                            |  
  ! |   we consider three different regions:                     |
  ! |    -  laurent extrapolated short range region (i)          |  
  ! |    -  spline-interpolated medium range region (ii)         |
  ! |    -  inverse-power extrapolated long range region (iii)   |
  ! |                                                            |
  ! |    matching of the first derivatives at the boundary       | 
  ! +------------------------------------------------------------+
  
  ! allocate the space for vlambda:
  allocate (v_lambda(num_pot_points, num_lambda_terms, &
       num_vib_levels, num_vib_levels) ) 
  allocate (spline_coeff(num_pot_points, num_lambda_terms, & 
       num_vib_levels, num_vib_levels) )
  allocate (r(num_pot_points))
  allocate (extrapolation_coeff(4, num_lambda_terms, & 
       num_vib_levels, num_vib_levels) )

  ! all processes allocate some scratch vectors needed in
  ! vpot(). These vectors are shared using potential_scratch
  ! module and are of constant size for a certain calulation
  allocate(v_tmp(num_pot_points))
  allocate(v_tmp2(num_pot_points))
  allocate(vtt(num_lambda_terms, num_vib_levels, num_vib_levels))

  ! reads the v_lambda(n,n') couplings from std input as
  ! generated by pot_gen.f program.
  rank_0_pot: if (myrank==0)then
     outer_vib_lev: do n=1,num_vib_levels
        inner_vib_lev: do nn=n,num_vib_levels
           read(potential_unit,*)i_dummy,i_dummy
           lambda_loop: do l=1,num_lambda_terms
              read(potential_unit,*)i_dummy
              potential_points_loop: do ir=1,num_pot_points
                 read(potential_unit,*,IOSTAT=ierr) r(ir), v
                 if(ierr>0) then
                    write(*,*)'something weird happened to input file' &
                         ,potential_unit
                    write(*,*)'program stopped'

#ifdef PARALLEL
                    write(*,*)'stopping mpi subsystem'
                    call mpi_finalize(ierr_mpi)
#endif

                    stop
                 endif
                 if(ierr<0) then
                    write(*,*)'unexpected end-of-file in input file' &
                         ,potential_unit
                    write(*,*)'program stopped'

#ifdef PARALLEL
                    write(*,*)'stopping mpi subsystem'
                    call mpi_finalize(ierr_mpi)
#endif

                    stop
                 endif
                 r(ir)=r(ir)/au_to_angst
                 v_lambda(ir,l,n,nn)=v*cm_to_au
                 v_lambda(ir,l,nn,n)=v*cm_to_au
              enddo potential_points_loop
           enddo lambda_loop
        enddo inner_vib_lev
     enddo outer_vib_lev
     
     if (points_to_discard<0) then 
        write(*,*) 'Warning: points_to_discard must be >= 0'
     endif
     write(*,*)' '
     write(*,1000)'+------------------------------+'
     write(*,1000)'|     potential processing     |'
     write(*,1000)'+------------------------------+'
     write(*,*)' '
     if (.not.use_variable_phase) then
        write(*,1017)'Starting point requested=',rbegin*au_to_angst,& 
             'matching=', rswitch*au_to_angst
     else
        write(*,1017)'Starting point requested=',rbegin*au_to_angst,& 
          'matching=', rend*au_to_angst
     endif
     write(*,1017)'first potential r value=',r(1)*au_to_angst,'last one=',&
          r(num_pot_points)*au_to_angst
     write(*,*)' '
     ! see which kind of extrapolation has been chosen
     ! and confirm it on std output
     if (rbegin<r(1).and.rend>r(num_pot_points))  choice=1
     if (rbegin<r(1).and.rend<=r(num_pot_points)) choice=2
     if (rbegin>=r(1).and.rend>r(num_pot_points)) choice=3

     select case (choice)
     case (1)
        write(*,999)'main: rbegin < r(1) and rend > r(num_pot_points)' 
        write(*,999)'Potential will be evaluate'
        write(*,1018)'with Laurent extrapolation for r <',r(1)*au_to_angst
        write(*,1018)'an with inverse power fitting for r >',    &
             r(num_pot_points)*au_to_angst
     case (2)
        write(*,999)'main: rbegin < r(1)' 
        write(*,999)'Potential will be evaluate'
        write(*,1018)'with a Laurent extrapolation for r <',r(1)*au_to_angst
     case (3)
        write(*,999)'main: rend > r(num_pot_points)' 
        write(*,999)'Potential will be evaluate'
        write(*,1018)'with an inverse power fitting for r >', &
             r(num_pot_points)*au_to_angst
     case default
        write(*,999)'main: rbegin > r(1) and rend < r(num_pot_points)'
        write(*,999)'standard fitting procedure with spline'
     end select

1017 format(a,t25,f15.4,5x,a,t50,f15.4)
1018 format(a,f6.2,5x,a,f6.2)
     write(*,1000)' '
     write(*,1000)'+-----------------------------+'
     write(*,1000)'| vlambda coefficients stored |'
     write(*,1000)'+-----------------------------+'

     ! Next subroutine generates an analitic extrapolation
     ! for the potential:
     ! for short range uses the form: a/R +b*R
     ! for long  range uses the form: 
     !           c*R**-n + d*R**-m
     !                
     ! Store the value of a(lambda,n,n'), b(lambda,n,n'),
     !   c(lambda,n,n') and d(lambda,n,n') in vector extra; 
     !   this array is used by vpot (see module potential).
     ! For values of r.gt.r(2) and r.lt.r(num_pot_points-1) generates the
     !   vspl vector which contains the spline coefficients that
     !   will be used by vpot through array vspl.
     call vfitt(num_vib_levels)
  endif rank_0_pot

#ifdef PARALLEL
  !     +--------------------------------------+
  !     | mpi: send potential to all processes |
  !     +--------------------------------------+
  dim1=num_pot_points*num_lambda_terms*num_vib_levels*num_vib_levels
  dim2=4*num_lambda_terms*num_vib_levels*num_vib_levels
  call mpi_bcast(v_lambda, dim1, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(spline_coeff, dim1, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(extrapolation_coeff, dim2, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD, ierr_mpi)
  call mpi_bcast(r, num_pot_points, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD, ierr_mpi)
#endif
  
  ! +-------------------------------------+
  ! |     main cycles of the program      |
  ! +-------------------------------------+
  ! cicle over the energies: this is the outermost loop
  if (myrank==0) then
     write(*,*)' '
     write(*,1000)'+---------------------------+'
     write(*,1000)'| main: loop over energies: |'
     write(*,1000)'+---------------------------+'
  endif
  
  rend_backup=rend
  rswitch_backup=rswitch
  loop_over_energy: do en_index=1,num_of_energies
     ! we set first_time here:
     ! We calculate in close_coupling() the step_size for the 
     !   log-der propagator only for the
     !   first energy (calculation is more stable), the first 
     !   partial wave and the first symmetry for CC and CS  
     if (en_index==1) first_time=.true.

     rend=rend_backup
     rswitch=rswitch_backup
     if(logarithmic_e) then
        collision_energy=first_energy*delta_e**(en_index-1)
     else
        collision_energy=first_energy+delta_e*real(en_index-1,kind=wp)
     endif

     ! total energy: 
     ! the next amount of energy is present 
     ! in the singlet, doublet, and triplet case:
     
     total_energy = vibr_energy(initial_vib_level) + &
     rot_constant(initial_vib_level)*real(init_N_state,kind=wp) &
                                    *(real(init_N_state,kind=wp)+uno) 

    ! we have to add a little more if the spin state
    ! is different of singlet:

      two_cases: if(doublet)then
      !doublet 

         if (real(init_N_state,kind=wp)<xinit_j_state)then
          ! F1 : j=N+1/2 N=j-1/2 
          total_energy=total_energy + &
          spin_rot*.5d0*real(init_N_state,kind=wp)
         elseif(real(init_N_state,kind=wp)>xinit_j_state) then
          ! F2 : j=N-1/2 N=j+1/2
          total_energy=total_energy - &
          spin_rot*.5d0*(real(init_N_state,kind=wp)+uno)
         endif

      elseif(triplet)then
      !triplet

        if (real(init_N_state,kind=wp)<xinit_j_state)then
         ! F1 : j=N+1 N=j-1
         total_energy=total_energy + &
         rot_constant(initial_vib_level)* & 
         (dos*real(init_N_state,kind=wp)+tres) - &
         spin_spin - &
         sqrt(rot_constant(initial_vib_level)**2 * &
             (dos*real(init_N_state,kind=wp)+tres)**2 + &
             spin_spin**2-dos*spin_spin*rot_constant(initial_vib_level))+&
         spin_rot*(real(init_N_state,kind=wp)+uno)       
        elseif(real(init_N_state,kind=wp)>xinit_j_state) then
         ! F3: j=N-1 N=j+1
         total_energy=total_energy - &
         rot_constant(initial_vib_level)* &
         (dos*real(init_N_state,kind=wp)-uno) - &
         spin_spin + &
         sqrt(rot_constant(initial_vib_level)**2 * &
             (dos*real(init_N_state,kind=wp)-uno)**2 + &
             spin_spin**2-dos*spin_spin*rot_constant(initial_vib_level))-&
        spin_rot*real(init_N_state,kind=wp)
         else
         ! F2 : j=N1 N=j1
         total_energy=total_energy
        endif

      endif two_cases

     internal_energy=total_energy
     total_energy=total_energy+collision_energy

     ! reduced energy
     reduced_e = total_energy * two_mu

     if (myrank==0) then
        write(*,*)' '
        write(*,2020)'INTERNAL ENERGY (cm-1):', &
                          internal_energy/cm_to_au
        write(*,1020)'NO:',en_index,'COLLISION ENERGY (cm-1):', &
                      collision_energy/cm_to_au
        write(*,1020)'NO:',en_index,'TOTAL ENERGY (cm-1):', &
                      total_energy/cm_to_au
     endif

1020 format(a,i4,5x,a,es9.2)
1019 format(a,i4,a)
2020 format(12x,a,es9.2)
 
     total_xs2=0.0_wp

#ifdef PARALLEL
     jtotx=0
     jtotm=0
     jtot_per_proc=0
     !   +-------------------------------------------------------+
     !   | mpi: only root process does jtot loop treatmente then | 
     !   |      it sends the results to other processes          |
     !   +-------------------------------------------------------+
     rank0_proc: if (myrank==0) then
        write(*,*)' '
        write(*,fmt='(t15,a)')'       * MPI summary *        '
        write(*,fmt='(t1,a)')'-------------------------------------'

        ! for nproc<=4 we increase jtotmax to a multiple of nproc
        if (nproc<=4) then          
           ! strange case with 1 partial wave
           if ((jtotmax - jtotmin + 1)==1) then
              write(*,*)'You have chosen to process only 1 angular momentum'
              write(*,*)'parallelization is useless with one partial wave'
              jtotmax=jtotmax+(nproc-1)
             write(*,*)'I will increase jtotmax to the number of processors' 
           else
              write(*,*)'You have chosen to run on less that 4 processors'
              write(*,*)'I will increase jtotmax to a multiple of nproc'
           endif
           ! normal case
           n_pw = (jtotmax - jtotmin + 1)
           if (mod(n_pw,nproc) == 0) then  ! multiple of nproc
              pw_per_proc=n_pw/nproc
           else
              jtot_temp=( (n_pw/nproc) + 1 ) * nproc
              n_pw = (jtot_temp )
              pw_per_proc=n_pw/nproc
              write(*,fmt='(a,i4)')'jtotmax increased to', &
                                    jtot_temp+jtotmin-1
           endif
           do n=1,nproc
              jtot_per_proc(n)=pw_per_proc
              jtotm(n)=jtotmin + (n-1)*pw_per_proc
              jtotx(n)=jtotm(n) + (pw_per_proc-1)
           enddo
        endif

        ! the following code will work better with nproc.gt.4
        ! we assign the work to all the procs distributing njtot/nproc on
        ! every proc, and then the rest to the first nrest procs.
        if (nproc>4) then           
           ! strange case with 1 partial wave
           if ((jtotmax - jtotmin + 1)==1) then
              write(*,*)' '
              write(*,*)'You have chosen to process only 1 angular momentum'
              write(*,*)'parallelization is useless with one partial wave'
              jtotmax=jtotmax+(nproc-1)
             write(*,*)'I will increase jtotmax to the number of processors' 
           endif
           ! normal case
           n_pw= (jtotmax - jtotmin + 1)
           pw_per_proc= n_pw/nproc
           pw_remaining=mod(n_pw,nproc)
           i=0
           do n=1,nproc
              if(n<=pw_remaining) then 
                 jtot_per_proc(n)=pw_per_proc+1
                 i=i+jtot_per_proc(n) 
                 jtotm(n)=jtotmin + (n-1)*jtot_per_proc(n)
              else
                 jtot_per_proc(n)=pw_per_proc
                 i=i+jtot_per_proc(n)
                 jtotm(n)=jtotmin + i-1
              endif
              jtotx(n)=jtotm(n) + ( jtot_per_proc(n)-1) 
           enddo
        endif
        ! writing out a good summary of what will be performed
        write(*,fmt='(t1,a)')'------------------------------------------'
        write(*,661)'process ','assigned pw','from','to'
        write(*,fmt='(t1,a)')'------------------------------------------'
        do n=1,nproc
           write(*,660)n-1,jtot_per_proc(n),jtotm(n),jtotx(n)
        enddo
        write(*,*)' '
     endif rank0_proc
660  format(t1,i4,t20,i4,t40,i4,t50,i4)
661  format(t1,a,t20,a,t40,a,t50,a)
     call mpi_bcast( jtotm, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
     call mpi_bcast( jtotx, nproc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
#endif

     if (myrank==0) then
        write(*,*)' '
        write(*,1000)'+---------------------------+'
        write(*,1000)'|  loop over partial waves  |'
        write(*,1000)'+---------------------------+'
     endif

#ifdef PARALLEL
     call mpi_barrier(MPI_COMM_WORLD, ierr_mpi)
     ! jtotstep not yet implemented for parallel calculations
     if (myrank==0) then
        write(*,*)''
        write(*,*)'whatever jtotstep has been chosen, &
                   jtotstep=1 will be used'
        write(*,*)'  for parallel calculations, sorry'
        write(*,*)''
     endif
     jtotstep=1
#endif

#ifndef PARALLEL
     jtotm(1)=jtotmin; jtotx(1)=jtotmax
#endif

     loop_over_partial_waves: do jtot=jtotm(myrank+1),jtotx(myrank+1), &
                                      jtotstep

     sigma_JJ2(jtot)=0.0_wp

        ! find Omega values
        partial_xs2=0.0_wp
           ! for CC there are only two symmetries allowed
           ! 1 is for r+l odd, 2 is for r+l even
           omega_min=1
           omega_max=2

        loop_over_parity: do omega = omega_min, omega_max
           ! count local channels |l,j,n,r> (or |j,n> for CS)
           ! num_local_ch counts the channels according to omega and jtot
           num_local_ch=0
           if (debug>=2)then
              write(unit_mpi,*)' '
              write(unit_mpi,*)' '
              write(unit_mpi,*)'----------- counting channels -----------'
              write(unit_mpi,1031)'omega','jtot','j','lmin','lmax','N'
           endif
           do i=1,num_ch  ! these are the asymptotic channels
                 ! for CC calculates the allowed l values
                 j=j_qn(i)
                 kr=r_qn(i) 
                 lmin=abs(jtot-j)
                 lmax=jtot+j
                 if(debug>=2)then
                    write(unit_mpi,1032)omega,jtot,j,lmin,lmax,kr
                 endif
                 do l=lmin,lmax
                    ! l+j odd - iparity=1 
                    ! l+kr-2 odd --> iparity=1 --> same l+kr
                    if ( (omega==1).and.(.not.parity(l+kr)) ) then
                       num_local_ch=num_local_ch+1
                       if (debug>=2)then 
                    write(unit_mpi,fmt='(t30,a)')'*counted channel in sym1'
                       endif
                    endif
                    ! l+kr even --> iparity=0
                    if ( (omega==2).and.(parity(l+kr)) ) then
                       num_local_ch=num_local_ch+1
                       if (debug>=2)then
                    write(unit_mpi,fmt='(t30,a)')'*counted channel in sym2'
                       endif
                    endif
                 enddo
           enddo

           ! Here we have a value of num_local_ch that is the dimension
           ! of the scattering involved
           if (num_local_ch==0)  then 
              ! this happens when there are no channels
              ! normally for jtot=0 omega=1
              write(*,*)''
              write(*,fmt='(a,i3,1x,a,2x,f4.1,5x,a,2x,i3,5x,a)') &
                    'RK-',myrank,'JTOT= ',real(jtot)-corr, &
                   'Symmetry=',omega, 'no available channels'
              cycle loop_over_parity
           else
              write(*,*)''
              write(*,fmt='(a,i3,1x,a,2x,f4.1,5x,a,2x,i3,5x,a)') &
                    'RK-',myrank,'JTOT= ',real(jtot)-corr, &
                   'Symmetry=',omega, 'propagating solution'           
              write(*,fmt='(a,i3,1x,2(a,2x,i3,5x))')'RK-',myrank,& 
                   'Number of channels for symmetry ',&
                   omega,'is ',num_local_ch
           endif
           symmetry_xs2=0.0_wp          
 
           ! This is for calculating the step_size 
           ! for each energy or only the for the first 
           ! Jtot value only.
           ! if(jtot.eq.jtotm(myrank+1)) first_time=.true.
           if (debug>=1)then
              write(unit_mpi,1000)'+-------------------------------+'
              write(unit_mpi,1000)'| main: calling close_coupling  |'
              write(unit_mpi,1000)'+-------------------------------+'
           endif

           call close_coupling(jtot, steps, num_ch, open_ch, & 
                num_local_ch, omega, & 
                symmetry_xs2, myrank, nproc, first_time, &
                multiplicity)

           first_time=.false.

           ! check the num_local_ch flag: if set to zero this means that
           ! close_coupling() found no open channels for the propagation 
           ! and returned with no cross section
           if (num_local_ch==0) then
              write(*,fmt='(a,i3,1x,2(a,2x,i3,5x),a)') &
                    'RK-',myrank,'JTOT= ',jtot, &
                   'Symmetry=',omega, 'no open channels'
              cycle loop_over_parity
           endif

           if (debug>=1)then
              write(unit_mpi,1000)'+--------------------------------+'
              write(unit_mpi,1000)'|  main: exiting close_coupling  |'
              write(unit_mpi,1000)'+--------------------------------+'
           endif

           ! here we have on exit the symmetry_xs2 
           ! memory segment that contains
           ! the partial cross section for this jtot and omega:
           ! accumulate partial cross section in partial_xs2
           !
           ! The array "symmetry_xs2(n1,n2,j1,j2,r1,r2)" is the
           ! opacity function P_J(jp vp <-- j v) of Pack,
           ! page 635, formula (24).
           !
           ! The array "partial_xs2(n1,n2,j1,j2,r1,r2)" is the
           ! total cross section for J, i.e., each of the 
           ! terms in the sum of the total cross section in
           ! the expression (23), page 635, of Pack:
           ! partial_xs2(n1,n2,j1,j2,r1,r2) =
           ! sigma_J2(n1,n2,j1,j2,r1,r2,jtot) =
           ! (pi/k_jv**2)*(2*J+1)*P_J(jp vp <-- j v) 

           do i=1,open_ch 
              n1=n_qn(i)+1
              j1=j_qn(i)+1
              r1=r_qn(i)+2
              s1=s_qn(i)
              do ii=1,open_ch
                 n2=n_qn(ii)+1
                 j2=j_qn(ii)+1
                 r2=r_qn(ii)+2
                 s2=s_qn(ii)

                 xtmp2=symmetry_xs2(n1,n2,j1,j2,s1,s2)* &
                    pi/(reduced_e-e_level(i)*two_mu)*real(2*jtot+1,kind=wp)

                    partial_xs2(n1,n2,j1,j2,s1,s2)= &
                    partial_xs2(n1,n2,j1,j2,s1,s2)+xtmp2
  
                    sigma_J2(n1,n2,j1,j2,s1,s2,jtot)=xtmp2

                    sigma_J2(n1,n2,j1,j2,s1,s2,jtot)= &
                    sigma_J2(n1,n2,j1,j2,s1,s2,jtot)*(au_to_angst**2) 

                    sigma_JJ2(jtot)=sigma_JJ2(jtot)+ &
                    sigma_J2(n1,n2,j1,j2,s1,s2,jtot)

              enddo
           enddo
        enddo loop_over_parity

1031    format(t10,a,t20,a,t30,a,t40,a,t50,a,t60,a)
1032    format(t10,i3,t20,i3,t30,i3,t40,i3,t50,i3,t60,i3)
        ! we have that the memory segment partial_xs contains the 
        ! partial cross section (already averaged over parity loop).
        ! we accumulate results in total_xs.

1021    format(2x,a,t7,a,t12,a,t16,a,t25,a,t30,a,t35,a,t39,a,t44,a,10x,a)
1022    format(f4.1,t7,i3,t12,a,t16,i3,t25,i3,t30, &
               f4.1,t35,a,t39,i3,t44,f4.1,10x,es11.4)

     write(*,*)''
     write(unit_mpi,1040)'Partial cross section for jtot:', real(jtot)-corr
     write(*,*)''
     write(unit_mpi,1021)'J',' v1','*->',' v2',' N1',' j1', &
                                  '*-> ',' N2',' j2','    sigma_J'
      do i=1,open_ch
      n1=n_qn(i)+1
      j1=j_qn(i)+1
      r1=r_qn(i)+2
      s1=s_qn(i)
       do ii=1,open_ch
       n2=n_qn(ii)+1
       j2=j_qn(ii)+1
       r2=r_qn(ii)+2
       s2=s_qn(ii)
        xtmp2=partial_xs2(n1,n2,j1,j2,s1,s2)*(au_to_angst**2)
         if (xtmp2.gt.0.0_wp)then
         write(unit_mpi,1022)real(jtot)-corr,n1-1,'*->',n2-1,r1-2, &
                       real(j1-1,kind=wp)-corr, &
                       '*->',r2-2,real(j2-1,kind=wp)-corr,xtmp2
         endif

      total_xs2(n1,n2,j1,j2,s1,s2)=total_xs2(n1,n2,j1,j2,s1,s2) &
                                  +xtmp2*real(jtotstep,kind=wp)
       enddo
      enddo
 
        if(debug>=1) then
           write(unit_mpi,*)'-------------------------------------------'
        endif 
     enddo loop_over_partial_waves

#ifdef PARALLEL
     call mpi_barrier(MPI_COMM_WORLD, ierr_mpi)
#endif
     write(*,*)'' 
     write(*,966)'(Process', myrank, 'of', nproc, 'finished pw_loop)'
966  format(A,1X,I2,1X,A,1X,I2,1X,A)
     if (myrank==0) then
        write(*,*)' '
        write(*,1000)'+-----------------------------+'
        write(*,1000)'|  main:  jtot loop finished  |'
        write(*,1000)'+-----------------------------+'
     endif

#ifdef PARALLEL
     ! we collect the result gained by different processes into root 
     !  process and we sum them up
     dim1=num_vib_levels**2*(num_rot_levels)**2*(3)**2
     call mpi_reduce(total_xs2, real_total2, dim1, MPI_DOUBLE_PRECISION, &
          MPI_SUM, 0, MPI_COMM_WORLD, ierr_mpi)

#endif

    if (myrank==0) then
       write(*,*)''
       write(*,*)'+-----------------------------------------------------+'
       write(*,*)'+-----------------------------------------------------+'
       write(*,*)''
    endif

     rank0_xss: if (myrank==0) then
        write(*,fmt='(a)')'Integral cross section:'
        write(*,*)''
        write(*,1024)'E_coll(cm-1)','v1','*->', &
                     ' v2',' N1',' j1','*->',' N2',' j2','  sigma'
      do i=1,open_ch
      n1=n_qn(i)+1
      j1=j_qn(i)+1
      r1=r_qn(i)+2
      s1=s_qn(i)
       do ii=1,open_ch
       n2=n_qn(ii)+1
       j2=j_qn(ii)+1
       r2=r_qn(ii)+2
       s2=s_qn(ii)

#ifdef PARALLEL
                    xtmp2=real_total2(n1,n2,j1,j2,s1,s2)
#endif

#ifndef PARALLEL
                    xtmp2=total_xs2(n1,n2,j1,j2,s1,s2)
#endif
                    if (xtmp2>0.0_wp)then
                       write(*,1025)collision_energy/cm_to_au, &
                       n1-1,'*->',n2-1,r1-2,real(j1-1,kind=wp)-corr, &
                       '*->',r2-2,real(j2-1,kind=wp)-corr,xtmp2
                    endif
      enddo
     enddo
     endif rank0_xss

     if (myrank==0) then
        write(*,*)''
        write(*,*)'+-----------------------------------------------------+'
        write(*,*)'+-----------------------------------------------------+'
        write(*,*)''
     endif

! write rotationally summed cross sections from j=xjout

   sigma_el2=0.0_wp
   sigma_in2=0.0_wp

     if (myrank==0) then
      xxtmp2=0.0_wp
      do i=1,open_ch
      n1=n_qn(i)+1
      j1=j_qn(i)+1
      r1=r_qn(i)+2
      s1=s_qn(i)
       do ii=1,open_ch
         n2=n_qn(ii)+1
         j2=j_qn(ii)+1
         r2=r_qn(ii)+2
         s2=s_qn(ii)

                 if(singlet) then ! sing-doub-trip

  jout_cy1: if(j1==1+nint(xjout))then

#ifdef PARALLEL
               xxtmp2(n1,n2)=xxtmp2(n1,n2)+real_total2(n1,n2,j1,j2,s1,s2)
               if(n1==n2) then
               sigma_el2(n1)=sigma_el2(n1)+real_total2(n1,n2,j1,j2,s1,s2)
               elseif(n1/=n2) then
               sigma_in2(n1)=sigma_in2(n1)+real_total2(n1,n2,j1,j2,s1,s2) 
               endif
#endif

#ifndef PARALLEL
               xxtmp2(n1,n2)=xxtmp2(n1,n2)+total_xs2(n1,n2,j1,j2,s1,s2)
               if(n1==n2) then
               sigma_el2(n1)=sigma_el2(n1)+total_xs2(n1,n2,j1,j2,s1,s2)
               elseif(n1/=n2) then
               sigma_in2(n1)=sigma_in2(n1)+total_xs2(n1,n2,j1,j2,s1,s2)
               endif
#endif

                 endif jout_cy1

              elseif(doublet) then ! sing-doub-trip

  jout_cy22: if(j1==1+nint(xjout))then

#ifdef PARALLEL
               xxtmp2(n1,n2)=xxtmp2(n1,n2)+real_total2(n1,n2,j1,j2,s1,s2)
               if(n1==n2) then
               sigma_el2(n1)=sigma_el2(n1)+real_total2(n1,n2,j1,j2,s1,s2)
               elseif(n1/=n2) then
               sigma_in2(n1)=sigma_in2(n1)+real_total2(n1,n2,j1,j2,s1,s2)
               endif
#endif

#ifndef PARALLEL
               xxtmp2(n1,n2)=xxtmp2(n1,n2)+total_xs2(n1,n2,j1,j2,s1,s2)
               if(n1==n2) then
               sigma_el2(n1)=sigma_el2(n1)+total_xs2(n1,n2,j1,j2,s1,s2)
               elseif(n1/=n2) then
               sigma_in2(n1)=sigma_in2(n1)+total_xs2(n1,n2,j1,j2,s1,s2)
               endif
#endif

                 endif jout_cy22

                 elseif(triplet) then ! sing-doub-trip

  jout_cy33: if(j1==1+nint(xjout))then

#ifdef PARALLEL
               xxtmp2(n1,n2)=xxtmp2(n1,n2)+real_total2(n1,n2,j1,j2,s1,s2)
               if(n1==n2) then
               sigma_el2(n1)=sigma_el2(n1)+real_total2(n1,n2,j1,j2,s1,s2)
               elseif(n1/=n2) then
               sigma_in2(n1)=sigma_in2(n1)+real_total2(n1,n2,j1,j2,s1,s2)
               endif
#endif

#ifndef PARALLEL
               xxtmp2(n1,n2)=xxtmp2(n1,n2)+total_xs2(n1,n2,j1,j2,s1,s2)
               if(n1==n2) then
               sigma_el2(n1)=sigma_el2(n1)+total_xs2(n1,n2,j1,j2,s1,s2)
               elseif(n1/=n2) then
               sigma_in2(n1)=sigma_in2(n1)+total_xs2(n1,n2,j1,j2,s1,s2)
               endif
#endif

            endif jout_cy33

              endif ! sing-doub-trip
           enddo
        enddo
     endif

    write(*,fmt='(a,f4.1)')'Rotationally summed xs from j= ', xjout
    write(*,*)''
    write(*,1026)'E_coll(cm-1)','v1','*->',' v2',' sigma_rot'
     do n1=1,initial_vib_level
      do n2=1,initial_vib_level
     write(*,1027)collision_energy/cm_to_au,n1-1,'*->',n2-1,xxtmp2(n1,n2)
      enddo
     enddo

    if (myrank==0) then
       write(*,*)''
       write(*,*)'+-----------------------------------------------------+'
       write(*,*)'+-----------------------------------------------------+'
       write(*,*)''
    endif

    write(*,*)'Elastic and inelastic total cross sections:'
    write(*,1035)'v','  sigma_el ','  sigma_in '
    do n=1,initial_vib_level
     write(*,1036)n-1,sigma_el2(n),sigma_in2(n)
    enddo

    if (myrank==0) then
       write(*,*)''
       write(*,*)'+-----------------------------------------------------+'
       write(*,*)'+-----------------------------------------------------+'
       write(*,*)''
    endif

    write(*,*)'Study for the sigma serie:'
    write(*,1037)'J','sigma_J'
    do jtot=jtotm(myrank+1),jtotx(myrank+1), jtotstep
     if(sigma_JJ2(jtot)>0.0_wp)then
     write(*,1029)real(jtot)-corr,sigma_JJ2(jtot)
     endif
    enddo

    if (myrank==0) then
       write(*,*)''
       write(*,*)'+-----------------------------------------------------+'
       write(*,*)'+-----------------------------------------------------+'
       write(*,*)''
    endif

  enddo loop_over_energy

1024 format(a,1x,t20,a,t25,a,t30,a,t35,a,t40,a,t45,a,t50,a,t55,a,5x,a)
1025 format(es11.2,t20,i3,t25,a,t30,i3,t35,i3,t40, &
            f4.1,t45,a,t50,i3,t55,f4.1,t62,es11.4)
1026 format(a,1x,t20,a,t25,a,t30,a,3x,a)
1027 format(es11.2,t20,i3,t25,a,t30,i3,t35,es11.4)
1029 format(9x,f4.1,25x,es11.4)
1035 format(t10,a,5x,a,5x,a)
1036 format(8x,i2,5x,es11.4,5x,es11.4) 
1037 format(11x,a,29x,a)
  if ((debug>=1).and.(myrank==0))then
     write(*,1000)'+------------------------------+'
     write(*,1000)'| main:  energy loop finished  |'
     write(*,1000)'+------------------------------+'
  endif

  ! summary for energies (rotationally summed cross section) goes here:

#ifdef PARALLEL
  call mpi_finalize(ierr_mpi)
#endif

  stop
end program scattering

