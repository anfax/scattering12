!!!*************************************************************
! 文件/File: eigenfull.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: eigenfull.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

subroutine eigenfull(Hmat,Omat,rowind,colptr,max_index,nnz,temp,eval,evec,nev,nconv,verbosity)
implicit none
  INTEGER, intent(IN) :: max_index,nnz,nev,verbosity
  INTEGER, intent(OUT) :: nconv
  INTEGER :: colptr(max_index+1), rowind(nnz),mtype,i,j
  REAL(kind=8) :: shift
  REAL(kind=8) :: Hmat(nnz), Omat(nnz)
  REAL(kind=8), TARGET :: eval(nev), evec(max_index,nev)

  REAL(kind=8) ::  dnrm2, dlapy2,temp
    external          dnrm2, dlapy2

  INTEGER :: iparam(11), ipntr(11)
  INTEGER :: ldv, ncv, info, ido, lworkl, maxitr, itercount, mxitercount, ierr,nrhs
  REAL(kind=8), ALLOCATABLE ::  a(:,:), b(:,:),w(:)

allocate(a(max_index,max_index))
allocate(b(max_index,max_index))
allocate(w(max_index))
a(:,:)=0.d0
b(:,:)=0.d0
do i=1,max_index
do j=colptr(i),colptr(i+1)-1
a(i,rowind(j))=Hmat(j)
b(i,rowind(j))=Omat(j)
end do
end do
!do i=1,max_index
!do j=1,max_index
!write(32,*)a(i,j),b(i,j)
!end do
!end do
!call flush_(32)

 call diagfull(a,b,max_index,max_index,w)
evec(:,:)=a(:,1:nev)
eval(:)=w(1:nev)
deallocate(a,b,w)
end subroutine eigenfull

 subroutine diagfull(a,b,n,lda,w)
    use nrtype, only : dbl,i4b
    use open_information
    implicit none
    integer(kind=i4b) itype,lda,ldb,info,lwork,i,j,ldvr,ldvl,nrhs,n,m,counter2
    character jobvl,jobvr,jobz,uplo
    real(kind=dbl),dimension(n,n) :: A
    real(kind=dbl),dimension(n,n) :: B
    real(kind=dbl),dimension(n) :: w
    real(kind=dbl),dimension(:,:),allocatable ::vl,vr
    real(kind=dbl),dimension(:),allocatable :: work,alphar,alphai,den,wr,wi
    counter2=counter2+1
    uplo='u'
    jobz='v'
    jobvl='n'
    jobvr='v'
    lda=n
    ldb=lda
    ldvr=n
    ldvl=n
    lwork=3*n-1
    itype=1
    m=lda
    ldb=lda
    allocate(work(lwork))
    call dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
90  format(4d14.7)
    write(23,*)w
!call flush_(23)
    write(6,*)'info=',info
    deallocate(work)
  end subroutine diagfull
