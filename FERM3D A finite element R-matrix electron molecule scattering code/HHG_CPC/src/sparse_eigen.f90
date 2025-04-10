!!!*************************************************************
! 文件/File: sparse_eigen.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: sparse_eigen.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


module sparse1
  use nrtype, only : dbl,i4b
  real(kind=dbl),dimension(:),allocatable :: Gamma0sparse,Olap_Sparse
  integer(kind=i4b),dimension(:),allocatable :: rowind
end module sparse1

subroutine sparse_eigen(Nclosed,Nopen)
  use nrtype, only : dbl,i4b
  use Brep,only : p_Patch,R0,x_coord=>coord_x,y_coord=>coord_y,z_coord=>coord_z
  use Open_information, only : bionda2_boundary,open_vect,max_index
  use control, only : option_solver,option_diagonalizer
  use Matrices, only : area_overlap_v=>area_overlap,area_overlap_c,gamma_c,gamma_v=>gamma,nnz
  use sparse1
  implicit none
  integer(kind=i4b) :: Row,Col,i,j,nj,dim,NOpen,Nclosed,info,LDB,counter,q
  integer(kind=i4b),parameter :: ni=216
  real(kind=dbl) :: tolerance,temp
  integer(kind=i4b),dimension(:),allocatable :: colptr,newrowind,newcolptr
  real(kind=dbl),dimension(:,:),allocatable :: olap
  real(kind=dbl),dimension(:),allocatable :: eval
  integer(kind=i4b) ::  nev,verbosity,nconv
  write(6,*)'bestia'
  nev=30
  temp= -50.0d0
  nconv=0
  verbosity=3
  counter=counter+1
  tolerance=1.0d-14
  nj=p_Patch
  LDB=Nclosed
  dim=ni*nclosed
  write(6,*)'before allocateion'
  allocate(Gamma0sparse(dim))
  allocate(Olap_Sparse(dim))
  allocate(rowind(dim+1))
  allocate(colptr(dim+1))
  allocate(eval(nev))
  !deallocate(open_vect)
  allocate(open_vect(NCLosed,nev))
  write(6,*)'before sort2'
  call sort2(p_Patch)
  write(6,*)'after sort2'
  !---------------------------
  do i=1,dim
     Gamma0sparse(i)=0.0d0
     Olap_sparse(i)=0.0d0
     rowind(i)=0
     colptr(i)=0
  end do
  Col=0 
  nnz=0
  !  colptr(0)=0
  write(6,*)'after initializing' 
  do j=1,nj
     if (bionda2_boundary(j).ne.'up_fun') then
        write(35,*)'up_fun'
        Col=Col+1
        colptr(Col+1)=0
        colptr(1) = 1
        Row=0
        do i=1,ni
           write(35,*)'up_fun2'
           if ((gamma_c(j,i).ne.0).and.(bionda2_boundary(gamma_c(j,i)).ne.'up_fun')) then
              Row=area_overlap_c(gamma_c(j,i))
              if ((dabs(gamma_v(j,i)).ge.tolerance).or.(dabs(area_overlap_v(j,i)).ge.tolerance)) then
                 write(35,*)'up_fun3'
                 nnz=nnz+1
                 !************************
                 Gamma0Sparse(nnz) = gamma_v(j,i)
                 !             Gamma0Sparse(nnz) = gamma_v(j,i)+area_overlap_v(j,i)
                 Olap_sparse(nnz) =  area_overlap_v(j,i)
                 !if (j.eq.gamma_c(j,i)) then
                 !             Olap_sparse(nnz) = 1.d0
                 !else 
                 !             Olap_sparse(nnz) = 0.d0
                 !end if
                 !************************
                 rowind(nnz)=Row
                 colptr(Col+1)=colptr(Col+1)+1
              end if
           end if
        enddo
        colptr(Col+1)=colptr(Col+1)+colptr(Col)
        !     write(60,*)colptr(Col),Col     
     end if
  enddo

910 format(f8.5)
  !write(362,*)Gamma0sparse
  !write(363,*)Olap_sparse
  !write(364,*)rowind
  !write(365,*)Colptr
  !if (counter.eq.1) then
  !write(469,*)RHS_co
  !write(470,*)rowind
  !write(471,*)Colptr
  !end iF
  write (6,*) "actual nnz=",nnz

  ! Find RHS_co = (Gamma_cc)**-1 * Gamma_co by solving the equivalent linear
  !  equation, Gamma_cc * RHS_co = Gamma_co
  !
  ! Note that only GammaSparse is packed in sparse format, and RHS_co is
  !  just a normal 2-d matrix, initially containing Gamma_co, that is
  !  overwritten with the solution
  !
  ! Note that LDB is in the row dimension of B (usually equal to nclosed)
  !
  ! dslusolve_(int n, int nnz, int nrhs, double values,
  !                int rowind, int colptr, double b, int ldb,
  !                int info, char save_factorization, int verbosity)

  !    call dslusolve(Nclosed,nnz,Nopen,Gamma_cc,rowind,colptr,RHS_co,LDB,info,'1',1)

  !write(366,*)Nclosed
  !write(366,*)nnz
  !write(366,*)Nopen
  !write(366,*)LDB
  !write(366,*)rowind
  !write(366,*)colptr
  !write(367,900)RHS_co
  !write(368,*)Gamma0Sparse
  write(6,*)'before eigensparse'
  deallocate(gamma_c)
  deallocate(area_overlap_v)
  deallocate(gamma_v)
!  call eigensparse(gamma0sparse,Olap_sparse,rowind,colptr,NClosed,nnz,temp,eval,open_vect,nev,nconv,verbosity)

  !call eigensparse_mein(colptr,NClosed,nnz,temp,eval,open_vect,nev,nconv,verbosity)
  !$!$   call dslusolve(Nclosed,nnz,Nopen,Gamma0sparse,rowind,colptr,RHS_co,LDB,info,'1',2)
  !  write(6,*)'mostro'
  !  if (counter.eq.2) then
!  if (option_solver.eq.'SuperLU') then
!     call eigensparse(gamma0sparse,Olap_sparse,rowind,colptr,max_index,nnz,temp,eval,open_vect,nev,nconv,verbosity)
!  else if ((option_solver.eq.'Pardiso').and.(option_diagonalizer.eq.'Arpack')) then
     ! 2) Option 2: use Pardiso for factorization
     ! Pardiso
!     call eigensparse_pardiso(gamma0sparse,Olap_sparse,rowind,colptr,max_index,nnz,temp,eval,open_vect,nev,nconv,verbosity)
!  else if ((option_solver.eq.'Pardiso').and.(option_diagonalizer.eq.'full')) then 
      call eigenfull(gamma0sparse,Olap_sparse,rowind,colptr,max_index,nnz,temp,eval,open_vect,nev,nconv,verbosity)
!  end if

900 format(f8.5)
  do j=1,nev
     q=0
     do i=1,NClosed,8
        q=q+1
        write(900+j,902)x_coord(q),y_coord(q),z_coord(q),open_vect(i,j)
     end do
  end do
  !  write(6,*)'info',verbosity
  write(38,*)'eval',eval
  !  end if
  deallocate(colptr)
902 format(4e13.5)
end subroutine sparse_eigen



