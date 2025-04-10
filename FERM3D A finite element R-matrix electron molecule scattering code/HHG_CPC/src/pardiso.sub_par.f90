!!!*************************************************************
! 文件/File: pardiso.sub_par.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: pardiso.sub_par.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

module pardiso_mod
use nrtype, only : i4b,dbl
contains
subroutine pardiso_sym(a,ja,ia,nnz,n,nrhs,rhs,ierr,pattern)
  IMPLICIT NONE
  INTEGER(kind=i4b) :: pt(64)
  integer(kind=i4b),intent(in) :: n,nnz,nrhs,pattern
  integer(kind=i4b) :: k,l,i,j,idum
  integer(kind=i4b),intent(inout):: ierr
  INTEGER(kind=i4b) :: maxfct, mnum, mtype, phase, error, msglvl
  INTEGER(kind=i4b) :: iparm(64)
  INTEGER(kind=i4b),intent(inout) :: IA(n+1), JA(nnz)
  real(kind=dbl),intent(inout) :: A(nnz),RHS(n,nrhs)
  real(kind=dbl),dimension(:,:),allocatable :: x
  REAL(kind=dbl) waltime1, waltime2, ddum
  DATA  maxfct /1/, mnum /1/
  integer(kind=i4b) omp_get_max_threads
  logical ierr_check, pattern_check
  external omp_get_max_threads
  if (pattern.ne.11)print *, 'enter pardiso subroutine'
  allocate(x(n,nrhs))
  !C..
  !C.. Set up PARDISO control parameter
  !C..
  !---------------------------------
  ! Modification of zero-based arrays
  !write(6,*)'inside pardiso',ja(1)
  if (ja(1).eq.0) then
     ja(:)=ja(:)+1
     ia(:)=ia(:)+1
  end if
  !---------------------------------
  x(:,:)=0.d0
  ! Diagnostics
  !write(400,*)a(:)
  !write(401,*)ia(:)
  !write(402,*)ja(:)
  !write(403,*)rhs(:,:)
  !write(404,*)x(:,:)
  !write(405,*)n,nnz,nrhs,omp_get_max_threads()
  !------------------------
  do i = 1, 64
     iparm(i) = 0
  end do
  iparm(1) = 1 ! no solver default
  iparm(2) = 2 ! fill-in reordering from METIS
  iparm(3) = omp_get_max_threads() ! numbers of processors, value of OMP_NUM_THREADS
  iparm(4) = 0 ! no iterative-direct algorithm
  iparm(5) = 0 ! no user fill-in reducing permutation
  iparm(6) = 0 ! =0 solution on the first n compoments of x
  iparm(7) = 0 ! not in use
  iparm(8) = 9 ! numbers of iterative refinement steps
  iparm(9) = 0 ! not in use
  iparm(10) = 13 ! perturbe the pivot elements with 1E-13
  iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
  iparm(12) = 0 ! not in use
  iparm(13) = 0 ! not in use
  iparm(14) = 0 ! Output: number of perturbed pivots
  iparm(15) = 0 ! not in use
  iparm(16) = 0 ! not in use
  iparm(17) = 0 ! not in use
  iparm(18) = -1 ! Output: number of nonzeros in the factor LU
  iparm(19) = -1 ! Output: Mflops for LU factorization
  iparm(20) = 0 ! Output: Numbers of CG Iterations
  error = 0 ! initialize error flag
  msglvl = 1 ! print statistical information
  !pattern_check=PRESENT(pattern)
  !if (pattern_check) then
  !   mtype = 11 ! nonsymmetric, indefinite
  !else
  !   mtype = -2 ! symmetric,indefinite
  !end if
  mtype=pattern
  !C.. Initiliaze the internal solver memory pointer. This is only
  !C necessary for the FIRST call of the PARDISO solver.
  do i = 1, 64
     pt(i) = 0
  end do
  !write(6,*)'1'
  !write(34,*)a
  !write(35,*)ia
  !write(36,*)ja
  !write(37,*)pt, maxfct, mnum, mtype, phase, n, idum, nrhs, iparm, msglvl, ddum, ddum, error
  !C.. Reordering and Symbolic Factorization, This step also allocates
  !C all memory that is necessary for the factorization
  phase = 11 ! only reordering and symbolic factorization
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,&
       & idum, nrhs, iparm, msglvl, ddum, ddum, error)
  !stop
  if (mtype.ne.11)WRITE(*,*) 'Reordering completed ... '
  IF (error .NE. 0) THEN
     WRITE(*,*) 'The following ERROR was detected: ', error
     STOP
  END IF
   if (mtype.ne.11)WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
   if (mtype.ne.11)WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
   if (mtype.ne.11)write(6,*)'2',mtype
  !C.. Factorization.
  phase = 22 ! only factorization
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,&
       & idum, nrhs, iparm, msglvl, ddum, ddum, error)
  if (mtype.ne.11)WRITE(*,*) 'Factorization completed ... '
  IF (error .NE. 0) THEN
     WRITE(*,*) 'The following ERROR was detected: ', error
     STOP
  ENDIF
  !C.. Back substitution and iterative refinement
  iparm(8) = 2 ! max numbers of iterative refinement steps
  phase = 33 ! only factorization
  if ((.not.pattern_check).and.(error.ne.0)) write(6,*)'3',error
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,&
       & idum, nrhs, iparm, msglvl, rhs, x, error)
   if (mtype.ne.11)WRITE(*,*) 'Solve completed ... '
   if (mtype.ne.11)WRITE(*,*) 'The solution of the system is '
  !DO i = 1, n
  !   WRITE(*,*) ' x(',i,') = ', x(i,:)
  !END DO
  !C.. Termination and release of memory
  phase = -1 ! release internal memory
  if ((.not.pattern_check).and.(error.ne.0))write(6,*)'4',error
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,&
       & idum, nrhs, iparm, msglvl, ddum, ddum, error)
  !write(51,*)x(:,1)
  if ((.not.pattern_check).and.(error.ne.0))write(6,*)'5',error
  rhs(:,:)=x(:,:)
  !ierr_check=PRESENT(ierr)
  !if (ierr_check)
  ierr=error
  deallocate(x)
END subroutine pardiso_sym

end module pardiso_mod
