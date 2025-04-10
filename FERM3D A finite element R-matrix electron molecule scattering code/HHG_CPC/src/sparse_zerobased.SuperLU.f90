!!!*************************************************************
! 文件/File: sparse_zerobased.SuperLU.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: sparse_zerobased.SuperLU.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


subroutine sparse(Nclosed,Nopen)
  use nrtype, only : dbl,i4b
  use Brep,only : p_Patch
  use Open_information, only : bionda2_boundary
  use Matrices, only : area_overlap_v=>area_overlap,RHS_co=>d,area_overlap_c,gamma_c,gamma_v=>gamma,nnz
  use control, only : iter_cont,option_Solver
  use pardiso_mod
  implicit none
!  interface 
!     subroutine pardiso_sym(a,ja,ia,nnz,n,nrhs,rhs,ierr,pattern)
!       integer*4,intent(in) :: n,nnz,nrhs
!       integer*4,intent(inout),optional :: ierr,pattern 
!       INTEGER*4,intent(inout) :: IA(n+1), JA(nnz)
!       real*8,intent(inout) :: A(nnz),RHS(n,nrhs)
!     end subroutine pardiso_sym
!  end interface

  integer(kind=i4b) :: Row,Col,i,j,nj,dim,NOpen,Nclosed,info,LDB,counter,mtype,ierr
  integer(kind=i4b),parameter :: ni=216
  real(kind=dbl) :: tolerance
  integer(kind=i4b),dimension(:),allocatable :: rowind,colptr,newrowind,newcolptr
  real(kind=dbl),dimension(:,:),allocatable :: olap
  real(kind=dbl),dimension(:),allocatable :: Gamma0sparse,newGammasparse
  real(kind=dbl),dimension(:),allocatable :: matlab_row,matlab_col,diag,somma
  !real*4 dtime, delta, tarray(2)
  !external dtime
  !delta = dtime ( tarray )
  !***************
  !---------------------------
  ! COMMENT
  ! bridge to SuperLU linear solver
  ! Sparse solver setup, store sparse matrices in format useful for SuperLU 2.0
  ! (nonsymmetric compressed row)
  !---------------------------

  write(6,*)'sparse_sparse subroutine'
  counter=counter+1
  tolerance=1.0d-14
  nj=p_Patch
  LDB=Nclosed
  dim=ni*nj
  allocate(Gamma0Sparse(nnz))
  allocate(rowind(nnz))
  allocate(colptr(0:p_Patch))
  call sort2(p_Patch)
  !---------------------------
  do i=1,nnz
     Gamma0sparse(i)=0
     rowind(i)=0
  end do
  do i=0,p_Patch
     colptr(i)=0
  end do
  Col=0 
  nnz=0
  colptr(0)=0
  if (option_solver.ne.'SuperLU') then
     do j=1,nj
        if (bionda2_boundary(j).ne.'up_fun') then
           Col=Col+1
           colptr(Col)=0
           Row=0
           do i=1,ni
              if ((gamma_c(j,i).ne.0).and.(bionda2_boundary(gamma_c(j,i)).ne.'up_fun')) then
                 Row=area_overlap_c(gamma_c(j,i))
                 if ((dabs(gamma_v(j,i)).ge.tolerance).or.(dabs(area_overlap_v(j,i)).ge.tolerance)) then
                    !if (Row.ge.Col) then ! COMMENTED BECAUSE SuperLU 2.0 DOES NOT HAVE SYMMETRIC MODE
                       nnz=nnz+1
                       Gamma0Sparse(nnz) = gamma_v(j,i)+area_overlap_v(j,i)
                       rowind(nnz)=Row-1
                       colptr(Col)=colptr(Col)+1
                    !end if
                 end if
              end if
           enddo
           colptr(Col)=colptr(Col)+colptr(Col-1)
        end if
     enddo
  else
     do j=1,nj
        if (bionda2_boundary(j).ne.'up_fun') then
           Col=Col+1
           colptr(Col)=0
           Row=0
           do i=1,ni
              if ((gamma_c(j,i).ne.0).and.(bionda2_boundary(gamma_c(j,i)).ne.'up_fun')) then
                 Row=area_overlap_c(gamma_c(j,i))
                 if ((dabs(gamma_v(j,i)).ge.tolerance).or.(dabs(area_overlap_v(j,i)).ge.tolerance)) then
                    nnz=nnz+1
                    Gamma0Sparse(nnz) = gamma_v(j,i)+area_overlap_v(j,i)
                    rowind(nnz)=Row-1
                    colptr(Col)=colptr(Col)+1
                 end if
              end if
           enddo
           colptr(Col)=colptr(Col)+colptr(Col-1)
        end if
     enddo
  end if
  deallocate(area_overlap_v)

  write(6,*)'before pardiso'
  if (option_solver.eq.'SuperLU') then
!     call dslusolve(Nclosed,nnz,Nopen,Gamma0sparse,rowind,colptr(0),RHS_co,LDB,info,'1',3)
  else
     mtype=-2
     ierr=0
     call pardiso_sym(Gamma0sparse,rowind,Colptr(0),nnz,nClosed,nOpen,RHS_co,ierr,mtype)
     write(6,*)'sparse error',ierr
  end if

  !delta = dtime ( tarray )
  !write(37,*)'delta end sparse=',delta
  deallocate(Gamma0sparse)
  deallocate(rowind)
  deallocate(colptr)

end subroutine sparse

subroutine sort2(p_Patch)
  use Matrices , only : open_mat_c => gamma_c,overlap_v => area_overlap,open_mat_v=>gamma
  use nrtype
  implicit none
  integer(kind=i4b) counter
  real(kind=dbl) :: prov,provola
  integer(kind=i4b) :: i,j,inc,l,max,s,p_Patch
  integer(kind=i4b),parameter :: n=216
  real(kind=dbl) :: V
  save counter
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
              if (counter.eq.1) then  
                 prov=open_mat_v(s,i)
              end if
              provola=overlap_v(s,i)
              j=i
              do
                 if (open_mat_c(s,j-inc) <= v) exit
                 open_mat_c(s,j)=open_mat_c(s,j-inc)
                 if (counter.eq.1) then
                    open_mat_v(s,j)=open_mat_v(s,j-inc)
                 end if
                 overlap_v(s,j)=overlap_v(s,j-inc)
                 j=j-inc
                 if(j <= inc) exit
              end do
              open_mat_c(s,j)=v
              if (counter.eq.1) then 
                 open_mat_v(s,j)=prov
              end if
              overlap_v(s,j)=provola
           end if
        end do
        if (inc <= 1) exit
     end do
  end do
  do j=1, max
  end do
end subroutine sort2


