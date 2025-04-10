!!!*************************************************************
! 文件/File: diago_sparse.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: diago_sparse.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


subroutine diag_sparse(p_Patch)
  use nrtype, only : dbl,i4b
  use Open_information
  use over
  use control, only : option_solver,option_diagonalizer
  implicit none
  integer(kind=i4b) nj,ni,p_Patch,i,j,nnz,dim,nev,nconv,m,l,Col,Row,verbosity,ind,row_sovrapp
  integer(kind=i4b),dimension(:),allocatable :: rowind,Colptr
  !-----
  !real(kind=dbl),dimension(max_index,max_index) ::open_mat_v,overlap_v 
  !-----
  real(kind=dbl),dimension(:),allocatable :: gamma0sparse,Olap_sparse,eval
  real(kind=dbl) temp
  real(kind=dbl),parameter :: tolerance =1.0d-14
  !real(kind=dbl),dimension(:,:),allocatable :: sovrapp_v

  ! To solve the eigenproblem: H*Psi = E*O*Psi
  !  where H is the Hamiltonian and O is the overlap matrix
  !
  ! call eigensparse( REAL H, REAL O, INT rowind, INT colptr,
  ! INT MatrixSize, INT nnz, REAL shift, REAL eval, REAL evec, INT nev,
  ! INT nconv, INT verbosity )
  !
  ! Note that H and O must be in arrays with identical sparsity structure
  !
  ! The variable "shift" is the initial guess for the iterative algorithm.
  ! Make sure this is just below the lowest eigenvalue you want to find.
  ! "nev" is the number of values that eigensparse will try to return.
  ! "nconv" is the number of returned values that are fully converged.
  ! (It's a good idea to compare nconv and nev afterward calling
  ! eigensparse!) "nnz" is the number of nonzeros, that is, the dimension of
  ! the sparse arrays. Note that rowind has dimension nnz, but colptr has
  ! dimension (0:MatrixSize).
  !
  ! Standard verbosity is 1. Setting verbosity to 2 or 3 gives additional
  ! information about normalization/convergence.
  !
  nev=min(max_index-1,(max_l+1)**2)  ! TEST
  verbosity=3
  nconv=0
  dim=max_index*36
  allocate(gamma0sparse(dim))
  allocate(Olap_sparse(dim))
  allocate(rowind(dim+1))
  allocate(ColPtr(dim+1))
  allocate(eval(nev))
  allocate(open_vect(max_index,nev))
  !***
  !allocate(sovrapp_v(max_index,max_index))
  !***
  !write(381,*)open_mat_v
  temp=-50.0d0
  nj=p_Patch
  ni=36
  do i=1,dim
     Gamma0sparse(i)=0.0d0
     rowind(i)=0
     colptr(i)=0
  end do
  do i=1,nev
     eval(i)=0.0d0
  end do
  do i=1,max_index
     do j=1,nev
        open_vect(i,j)=0.0d0
     end do
  end do
  !--------------------------
  ! Compressed row format
  !--------------------------
  nnz=0
  !  colptr(0)=0
  Col=0
  ind=0
  do j=1,nj
     if (bionda2_boundary(j).eq.'up_fun') then
        ind=ind+1
        overlap_c(j,1)=ind
     end if
  end do

  do j=1,nj
     if (bionda2_boundary(j).eq.'up_fun') then
        Col=Col+1
        colptr(Col+1)=0
        colptr(1) = 1
        Row=0
        row_sovrapp=0
        do i=1,ni
           if ((open_mat_c(j,i).ne.0).and.(bionda2_boundary(open_mat_c(j,i)).eq.'up_fun')) then
              Row=overlap_c(open_mat_c(j,i),1)
              nnz=nnz+1
              Gamma0Sparse(nnz) = open_mat_v(j,i)
              Olap_sparse(nnz) = overlap_v(j,i)
              !***
              !sovrapp_v(Col,Row) = overlap_v(j,i)
              !***
              rowind(nnz)=Row
              colptr(Col+1)=colptr(Col+1)+1
           end if
        enddo
        colptr(Col+1)=colptr(Col+1)+colptr(Col)
     end if
  enddo
  !$write(18,*)'nnz',nnz
  deallocate(open_mat_c)
  deallocate(open_mat_v)
  deallocate(overlap_v)
  !do i=1,nnz
  !write(3336,9000)Gamma0sparse(i),Olap_sparse(i),rowind(i),i
  !end do
9000 format(2f14.7,2I6)
  !--------------------------------------------------------------------------------------------
  ! Solution of the eigenvalue problem
  ! 1) Option 1: use SuperLU for factorization
  ! SuperLU
!  if (option_solver.eq.'SuperLU') then
!     call eigensparse(gamma0sparse,Olap_sparse,rowind,colptr,max_index,nnz,temp,eval,open_vect,nev,nconv,verbosity)
!  else if ((option_solver.eq.'Pardiso').and.(option_diagonalizer.eq.'Arpack')) then
     ! 2) Option 2: use Pardiso for factorization
     ! Pardiso
!     call eigensparse_pardiso(gamma0sparse,Olap_sparse,rowind,colptr,max_index,nnz,temp,eval,open_vect,nev,nconv,verbosity)
!  else if ((option_solver.eq.'Pardiso').and.(option_diagonalizer.eq.'full')) then 
      call eigenfull(gamma0sparse,Olap_sparse,rowind,colptr,max_index,nnz,temp,eval,open_vect,nev,nconv,verbosity)
!  end if
  !--------------------------------------------------------------------------------------------
  deallocate(gamma0sparse)
  deallocate(Olap_sparse)
  deallocate(rowind)
  deallocate(colptr)
  !write(336,*)eval
  deallocate(eval)

  if (option_diagonalizer.ne.'full') write(6,*)verbosity

  !do j=1,nev
  !do i=1,max_index,4
  !  write(800+j,*)open_vect(i,j)
  !end do
  !end do
  !***
  !call order(open_vect)
  !***
  !$write(339,*)max_index
  !$write(338,*)open_vect
end subroutine diag_sparse


subroutine order(z)
  use open_information
  !use Solve
  use over
  use nrtype, only : dbl,i4b
  implicit none
  integer(kind=i4b) i,j,l,f,n
  real(kind=dbl),dimension(:,:),allocatable :: vl,vr,orto
  real(kind=dbl),dimension(max_index,max_index) :: z,finto
  !-------------------------------------------------------
  !Calculation of the scalar product of the open functions
  !-------------------------------------------------------
  n=max_index-1
  finto=sovrapp_v
  allocate(orto(max_index,max_index))
  do l=1,n
     do f=1,n
        orto(l,f)=0.0d0
        do i=1,n+1
           do j=1,n+1
              !write(201,*)orto(l,f),z(i,l),finto(i,j),z(j,f)
              orto(l,f)=orto(l,f)+z(i,l)*finto(i,j)*z(j,f)
           end do
        end do
        if (orto(l,f).lt. 0.0d0) then
           orto(l,f) = -orto(l,f)
        end if
        write(371,*)orto(l,f),l,f
     end do
  end do
  !write(376,*)z
  !--------------------------------------------------------
  !Orthogonalization and sorting of the functions
  !--------------------------------------------------------
  call ortog(z,finto,max_index,orto)
  !write(377,*)z
  !--------------------------------------------------------
  !Control of the orthogonalization
  !--------------------------------------------------------
  orto=0.0d0
  do l=1,n
     do f=1,n
        do i=1,n
           do j=1,n
              orto(l,f)=orto(l,f)+z(i,l)*finto(i,j)*z(j,f)
           end do
        end do
        if (orto(l,f).lt. 0.0d0) then
           orto(l,f) = - orto(l,f)
        end if
        !4             write(373,*)orto(l,f),l,f
     end do
  end do
  allocate(open_vect(n,n))
  open_vect=z
  !      write(375,*)open_vect
end subroutine order



subroutine sort3(p_Patch)
  use Matrices , only : open_mat_c => gamma_c
  use nrtype
  implicit none
  integer(kind=i4b) counter
  real(kind=dbl) :: prov,provola
  integer(kind=i4b) :: i,j,inc,l,max,s,p_Patch
  integer(kind=i4b),parameter :: n=216
  real(kind=dbl) :: V
  save counter
  write(6,*)'sort3'
  counter=counter+1
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
              j=i
              do
                 if (open_mat_c(s,j-inc) <= v) exit
                 open_mat_c(s,j)=open_mat_c(s,j-inc)
                 j=j-inc
                 if(j <= inc) exit
              end do
              open_mat_c(s,j)=v
           end if
        end do
        if (inc <= 1) exit
     end do
  end do
  do j=1, max
  end do
end subroutine sort3

