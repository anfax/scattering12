!!!*************************************************************
! 文件/File: solve.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: solve.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


module Solve
  use nrtype, only : dbl,i4b
  real(kind=dbl),dimension(:),allocatable :: eigen_value
  real(kind=dbl),dimension(:,:),allocatable :: eigen_vector 
  real(kind=dbl) :: Z
  interface lin
     module procedure lin1,lin2,lin4
  end interface
  integer(kind=i4b) :: rowA

  ! interface invert
  !   module procedure invert_real,invert_complex
  ! end interface invert

  !----------------------------
  ! This module contains all the subroutines that perform the linear algebra
  ! tasks or setup them
  !----------------------------

contains
  subroutine lin1(gamma)
    use nrtype, only : dbl,i4b
    use Brep
    use control  
    implicit none
    real(kind=dbl),dimension(:,:),allocatable :: A,c
    real(kind=dbl),dimension(:),allocatable :: B,gamma_cond
    real(kind=dbl),dimension(p_Patch,p_Patch) :: gamma
    integer(kind=i4b) :: counter2,columnA,columnB,i,j,globnode,node
    real(kind=dbl),dimension(:),allocatable :: bound_cond
    allocate(gamma_cond(p_Patch))
    allocate(A(p_Patch,p_Patch))
    allocate(B(p_Patch))
    allocate(bound_cond(p_Patch)) 
    !---------------------------------------------------------------------
    !This subroutine solves the Laplace equation (on a rectangular membrane)
    !---------------------------------------------------------------------
    rowA=0
    columnA=0
    columnB=0
    counter2=0 
    ! specify the kind of problem that the program solves
    ! boundary value problem
    !------------------------  
    !set boundary conditions
    do node=1,p_Patch
       if (edge(node).eq.'up') then
          bound_cond(node)=0
       else if (edge(node).eq.'down') then
          bound_cond(node)=200
       else if (edge(node).eq.'left') then
          bound_cond(node)=0
       else if (edge(node).eq.'right') then
          bound_cond(node)=0
       end if
    end do
    !-------------------------- 
    do i=1,p_Patch
       do j=1,p_Patch
          if (i.eq.j) then
             !             write(128,*)gamma(i,j)
          end if
       end do
    enddo
    do node=1,p_Patch
       gamma_cond(node)=0.0d0
       B(node)=0.0d0
       do j=1,p_Patch
          A(node,j)=0.0d0
          if((edge(node).eq.'bulk')) then
             gamma_cond(node)=gamma_cond(node)-gamma(node,j)*bound_cond(j)
          end if
       end do
    end do
    do node=1,p_Patch
       counter2=counter2+1
       !       write(6,*)'counter2=',counter2
       if (edge(node).eq.'bulk') then
          rowA=rowA+1
          columnA=0
          do i=1,p_Patch
             if (edge(i).eq.'bulk') then
                columnA=columnA+1
                A(rowA,columnA)=gamma(node,i)
             end if
          end do
       end if
    end do
    i=0
    do node=1,p_Patch
       if (edge(node).eq.'bulk') then
          i=1+i
          if (i.eq.(rowA+1)) exit
          B(i)=B(i)+gamma_cond(node)
       end if
    end do
    deallocate(bound_cond)
    allocate(c(p_Patch,1))
    !    write(6,*)'donne'
    !    write(29,*)'b=',b
    do i=1,rowA
       c(i,1)=b(i)
    end do

    do i=1,p_Patch
       do j=1,p_Patch
          if (i.eq.j) then
             !             write(128,*)a(i,j)
          end if
       end do
    enddo
    call lin_solve(a,c,p_patch,p_patch,1)  
    !    deallocate(a)
    deallocate(c)
  end subroutine lin1



  subroutine lin2
    use nrtype, only: dbl,i4b
    use Brep
    use control  
    use Matrices
    use Open_information
    implicit none
    real(kind=dbl),dimension(:,:),allocatable :: c,parz
    real(kind=dbl),dimension(:),allocatable :: JJ,JJp,NN,NNp

    real(kind=dbl),dimension(max_open,max_open) :: surface_open_overlap,omega
    integer(kind=i4b) :: columnA,i,j,globnode,node,m,mm,mostro
    real(kind=dbl) :: sparseness_rate,zero
    !    allocate(A(p_Patch,64+max_open))
    rowA=0 
    columnA=0     
    !---------------------------------------------------------------------
    ! COMMENT
    ! This subroutine performs the linear algebra for the R-matrix, inversion of the closed
    ! closed part and diagonalization of the open-open part, as in Greene et
    ! al. Rev. Mod. Phys. 1996 
    !---------------------------------------------------------------------
    !XXXXXXXXXXXXXX
    if (calculation_type.eq.'bound') then
       ! Eigenvalue problem solution for bound states
       call sparse_eigen(p_Patch-max_index,max_open)
       stop 
    end if
    !XXXXXXXXXXXXXXXX
    !write(6,*)'rowA',rowA
    !    pause
    zero=0
    do i=1,rowA
       do j=1,rowA
          !          if (A(i,j).lt.(1.0d-8)) then
          !             zero=zero+1
          !          end if
       end do
    enddo
    !write(6,*)'zero',zero 
    !    sparseness_rate=zero/(rowA*rowA)
    !    write(6,*)'sparseness=',sparseness_rate 
    !    write(6,*)'giorgio2',rowA
    !------------------------------------------------------------------------
    !Introduce the matrix of the overlap between the harmonics on the surface
    !It's the identity  matrix, dimension (max_open,max_open)
    !------------------------------------------------------------------------
    do m=1,max_open
       do mm=1,max_open
          if (m.eq.mm) then
             if (m.eq.1) then
                surface_open_overlap(m,mm) = 1
             else 
                surface_open_overlap(m,mm) = 1
             end if
          else
             surface_open_overlap(m,mm) = 0
          end if
       end do
    end do

    !------------------------------------------------------------------------
    mostro=2
    select case(mostro)
    case(1)

       !------------------------------------------------------------------------
       !  Diagonalize omega (as in orange review)
       !------------------------------------------------------------------------


    case(2)
       !------------------------------------------------------------------------
       ! Use sparse matrix solver SuperLU or Pardiso
       !------------------------------------------------------------------------

       allocate(d(p_Patch-max_index,max_open))
       allocate(c(max_open,p_Patch-max_index))
       d=0.0d0
       c=0.0d0 
       do m=1,max_open
          columnA=0
          do node=1,p_Patch
             if (bionda2_boundary(node).ne.'up_fun') then
                columnA=columnA+1
                d(columnA,m)=open_closed(node,m)+area_overlap_open_closed(node,m)
                c(m,columnA)=open_closed(node,m)+area_overlap_open_closed(node,m) 
                !write(3001,*)d(columnA,m)
                !write(3002,*)open_closed(node,m)
                !write(3003,*)area_overlap_open_closed(node,m)
             end if
          end do
       end do

       rowA=p_Patch-max_index
       !------------------------------------------
       ! Sets up the matrices in compressed column format for sparse solvers
       call sparse(rowA,max_open)
       !write(6,*)'after sparse'
       !-------------------------------------------
900    format(f8.5)
       allocate(parz(max_open,max_open))
       parz=0.0d0
       omega=0.0d0
       ! Setting up the reduced omega matrix
       do m=1,max_open
          do mm=1,max_open
             do i=1,rowA
                parz(m,mm)=parz(m,mm)+c(m,i)*d(i,mm)
             end do
             omega(m,mm)=open_open(m,mm)+area_overlap_open_open(m,mm)-parz(m,mm) 
          end do
       end do
       deallocate(area_overlap_open_open)
       deallocate(area_overlap_open_closed)
       deallocate(c)
       deallocate(parz)
       !       call diagol(omega,surface_open_overlap,max_open,max_open)
       allocate(eigen_vector(max_open,max_open))
       allocate(eigen_value(max_open))    
       ! Diagonalization of the omega matrix
       call diag(omega,surface_open_overlap,max_open,max_open) 
       allocate(JJ(max_open))
       allocate(JJp(max_open))
       allocate(NN(max_open))
       allocate(NNp(max_open))
      
       ! Call routine that calculates matching functions at R0
       call cross_section(JJ,JJp,NN,NNp)
       ! Call routine that calculates scattering observables
       call kmatrix(JJ,JJp,NN,NNp,Z,E)

       deallocate(JJ)
       deallocate(JJp)
       deallocate(NN)
       deallocate(NNp)

       deallocate(eigen_vector)
       deallocate(eigen_value)

       deallocate(d)
       write(6,*)'end loop for energy ',E

       !------------------------------------------------------------------------
    end select
  end subroutine lin2



  subroutine lin4(gamma,area_overlap,bionda_boundary)
    use nrtype, only: dbl,i4b
    use Brep
    use control  
    implicit none
    real(kind=dbl),dimension(:,:),allocatable :: A,B,c,d
    real(kind=dbl),dimension(:),allocatable :: gamma_cond
    real(kind=dbl),dimension(p_Patch,p_Patch) :: gamma,area_overlap
    integer(kind=i4b) :: counter2,columnA,columnB,i,j,globnode,node
    real(kind=dbl),dimension(:),allocatable :: bound_cond,overlap
    character*10,dimension(p_Patch) :: bionda_boundary   
    allocate(gamma_cond(p_Patch))
    allocate(A(p_Patch,p_Patch))
    allocate(B(p_Patch,p_patch))
    allocate(bound_cond(p_Patch))
    rowA=0 
    columnA=0
    columnB=0
    counter2=0 
    A=0
    B=0  
    !---------------------------------------------------------------------
    !This subroutine solves the motion on the circular boundary
    ! (and finds the surface harmonics)
    !---------------------------------------------------------------------
    do node=1,p_Patch
       counter2=counter2+1
       if (bionda_boundary(node).eq.'up_fun') then
          rowA=rowA+1
          columnA=0
          do i=1,p_Patch
             if (bionda_boundary(i).eq.'up_fun') then
                columnA=columnA+1
                A(rowA,columnA)=gamma(node,i)
                B(rowA,columnA)=area_overlap(node,i)
             end if
          end do
       end if
    end do

    deallocate(bound_cond)
    allocate(c(rowA,rowA))
    allocate(d(rowA,rowA))

    do i=1,rowA
       do j=1,rowA
          d(i,j)=b(i,j)
          c(i,j)=a(i,j)
       end do
    enddo
    call diagol(c,d,rowA,rowA)
    allocate(overlap(rowA)) 
    do i=1,rowA
       do j=1,rowA
          overlap(i)=c(i,j)*b(j,i)            
       end do
    end do
    deallocate(a)
    deallocate(b)
    deallocate(c)
    deallocate(d)
  end subroutine lin4




  subroutine lin_solve(a,b,n,lda,nrhs)
    use nrtype, only : dbl,i4b
    use control
    implicit none
    integer(kind=i4b) lda,ldb,info,i,j,nrhs,n,m
    character trans
    real(kind=dbl),dimension(n,n) :: A
    real(kind=dbl),dimension(n,1) :: B
    integer(kind=i4b),dimension(:),allocatable :: ipiv 
    !------------------------------------------------------------------------
    ! Subroutine used to solve the linear system that rises in Laplace/Poisson equation
    ! Electrostatic potential problems
    !------------------------------------------------------------------------

    lda=n
    trans='n'
    m=lda
    ldb=lda 
    allocate(ipiv(n)) 
    call dgetrf(m,n,a,lda,ipiv,info)
    call dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
90  format(4d14.7)
    deallocate(ipiv)
  end subroutine lin_solve


  subroutine invert_real(a,n,lda)
    use control
    implicit none
    integer(kind=i4b) lda,ldb,info,i,j,n,m,lwork
    real(kind=dbl),dimension(n,n) :: A
    integer(kind=i4b),dimension(:),allocatable :: ipiv 
    real(kind=dbl),dimension(:),allocatable :: work

    !------------------------------------------------------------------------
    ! Inverts the matrix omega
    !------------------------------------------------------------------------
    ! The dimensioning of lwork is not the best possible
    !------------------------------------------------------------------------
    lda=n
    m=lda
    ldb=lda 
    lwork=4*n
    allocate(ipiv(n)) 
    allocate(work(lwork)) 
    !write(6,*)m,n,lda
    call dgetrf(m,n,a,lda,ipiv,info)
    !    write(6,*)'after dgetrf'
    if (info.ne.0) write(6,*)'info invert_real 1=',info
    call dgetri(n,a,lda,ipiv,work,lwork,info)
    if (info.ne.0) write(6,*)'info invert_real 2=',info

    !    write(6,*)'after dgetrs  info=',info
90  format(4d14.7)
    deallocate(ipiv)
    deallocate(work)
    !write(6,*)'after invert real'
  end subroutine invert_real



  subroutine invert_complex(a,n,lda)
    use nrtype,only : dbl,i4b,dpc
    use control
    implicit none
    integer(kind=i4b) lda,ldb,info,i,j,n,m,lwork
    complex(kind=dpc),dimension(n,n) :: A
    integer(kind=i4b),dimension(:),allocatable :: ipiv
    complex(kind=dpc),dimension(:),allocatable :: work
    !------------------------------------------------------------------------
    ! Inverts the matrix omega
    !------------------------------------------------------------------------
    ! The dimensioning of lwork is not the best possible
    !------------------------------------------------------------------------
    lda=n
    m=lda
    ldb=lda
    lwork=4*n
    allocate(ipiv(n))
    allocate(work(lwork))
    call zgetrf(m,n,a,lda,ipiv,info)
    if (info.ne.0)write(6,*)'info invert complex 1 =',info
    call zgetri(n,a,lda,ipiv,work,lwork,info)
    if (info.ne.0)write(6,*)'info invert complex 2 =',info
90  format(4d14.7)
    deallocate(ipiv)
    deallocate(work)
  end subroutine invert_complex





  subroutine diag(a,b,n,lda)
    use nrtype, only : dbl,i4b
    !    use control
    use open_information 
    implicit none
    integer(kind=i4b) itype,lda,ldb,info,lwork,i,j,ldvr,ldvl,nrhs,n,m,counter2
    character jobvl,jobvr,jobz,uplo
    real(kind=dbl),dimension(n,n) :: A
    real(kind=dbl),dimension(n,n) :: B
    real(kind=dbl),dimension(:,:),allocatable ::vl,vr
    real(kind=dbl),dimension(:),allocatable :: work,alphar,alphai,den,wr,wi,w
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
    allocate(w(n)) 
    allocate(work(lwork))
    !    allocate(vr(n,n))
    !    allocate(vl(n,n))
    !    allocate(alphai(n))
    !    allocate(alphar(n))
    !    allocate(den(n))
    !    allocate(wr(n))
    !    allocate(wi(n)) 
    !    call dgegv(JOBVL,JOBVR,N,A,LDA,B,LDB,ALPHAR,ALPHAI,DEN,VL,LDVL,VR,LDVR,WORK,LWORK,INFO)
    !    call   DGEEV( JOBVL, JOBVR, N, b, LDA, WR, WI, VL, LDVL, VR, LDVR,WORK, LWORK, INFO )
    call dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)

!    allocate(eigen_vector(n,n))
!    allocate(eigen_value(n))

    eigen_vector=a
    eigen_value=w
90  format(4d14.7)
    !    write(34,*)a 
    write(23,*)w
    !    write(233,*)a
    if (info.ne.0)write(6,*)'info diag=',info
!    call cross_section
    deallocate(work)
    deallocate(w)
    !    deallocate(vr)
    !    deallocate(vl)
    !    deallocate(alphai)
    !    deallocate(alphar)
    !    deallocate(den)
!    deallocate(eigen_vector)
!    deallocate(eigen_value)
  end subroutine diag


  subroutine diag_complex(a,b,n,lda,eigen)
    use nrtype, only : dbl,i4b,dpc
    use control
    use open_information
    implicit none
    integer(kind=i4b) itype,lda,ldb,info,lwork,i,j,ldvr,ldvl,nrhs,n,m,counter2
    character jobvl,jobvr,jobz,uplo
    complex(kind=dpc),dimension(n,n) :: A
    complex(kind=dpc),dimension(n,n) :: B
    complex(kind=dpc),dimension(:,:),allocatable ::vl,vr
    complex(kind=dpc),dimension(:),allocatable :: work,wr,wi,w
    real(kind=dbl),dimension(:),allocatable :: alphar,alphai,den,rwork
    complex(kind=dpc),dimension(n) :: eigen 
    counter2=counter2+1
    !    counter2=2 
    uplo='u'
    jobz='v'
    jobvl='n'
    jobvr='v'
    lda=n
    ldb=lda
    ldvr=n
    ldvl=n
    lwork=2*n
    itype=1
    m=lda
    ldb=lda
    allocate(w(n))
    allocate(work(lwork))
    allocate(vr(n,n))
    allocate(vl(n,n))
    allocate(rwork(2*n))
    !allocate(alphai(n))
    !allocate(alphar(n))
    !allocate(den(n))
    !allocate(wr(n))
    !allocate(wi(n))
    !allocate(rwork(lwork))

    call ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK,LWORK, RWORK, INFO )

90  format(4d14.7)
    !    write(34,*)a
    !    write(35,*)w
    if (info.ne.0)write(6,*)'info diag_complex=',info
    eigen=w 
    deallocate(rwork)
    deallocate(work)
    deallocate(vr)
    deallocate(vl)
    deallocate(w)
    !deallocate(alphai)
    !deallocate(alphar)
    !deallocate(den)
  end subroutine diag_complex


  subroutine diag_symm(a,n,lda,eigen)
    use nrtype, only : dbl,i4b
    !    use control
    use open_information
    implicit none
    integer(kind=i4b) itype,lda,ldb,info,lwork,i,j,ldvr,ldvl,nrhs,n,m,counter2
    character jobvl,jobvr,jobz,uplo
    real(kind=dbl),dimension(n,n) :: A
    real(kind=dbl),dimension(:,:),allocatable ::vl,vr
    real(kind=dbl),dimension(:),allocatable :: work,alphar,alphai,den,wr,wi,w
    real(kind=dbl),dimension(n) :: eigen
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
    allocate(w(n))
    allocate(work(lwork))
    !    allocate(vr(n,n))
    !    allocate(vl(n,n))
    !    allocate(alphai(n))
    !    allocate(alphar(n))
    !    allocate(den(n))
    !    allocate(wr(n))
    !    allocate(wi(n))
    !write(370,*)a
    call dsyev(jobz,uplo,n,a,lda,w,work,lwork,info) 
    !    allocate(eigen_vector(n,n))
    !    allocate(eigen_value(n))

    !    eigen_vector=a
    eigen=w
90  format(4d14.7)
    !    write(34,*)a
    if (info.ne.0)write(6,*)'info diag_symm=',info
    !    call cross_section
    deallocate(work)
    deallocate(w)
    !    deallocate(vr)
    !    deallocate(vl)
    !    deallocate(alphai)
    !    deallocate(alphar)
    !    deallocate(den)
    !    deallocate(eigen_vector)
    !    deallocate(eigen_value)

  end subroutine diag_symm


  subroutine diagol(a,b,n,lda)
    use nrtype, only : dbl,i4b
    use open_information, only : open_vect,max_open,max_index,E
    implicit none
    integer(kind=i4b) itype,n,lda,ldb,info,lwork,i,j,l,matz,nm,counter2,f
    character uplo,jobz
    real(kind=dbl),dimension(n,n) :: a,b
    real(kind=dbl),dimension(:,:),allocatable :: vl,vr,z,finto,orto
    real(kind=dbl),dimension(:),allocatable :: w,work,alphar,alphai,den,eigen
    !    counter2=counter2+1
    counter2=2
    matz=1
    jobz='n'
    nm=n
    uplo='u'
    lda=n
    ldb=lda
    lwork=9*n+1
    itype=1
    allocate(eigen(n))
    allocate(work(n))
    allocate(w(n))   
    allocate(z(n,n))    
    allocate(vr(n,n))
    allocate(vl(n,n))
    allocate(alphai(n))
    allocate(alphar(n))
    allocate(den(n)) 
    write(6,*)'n=',n
    !write(6,*)'antonio' 
    !    if (counter2.eq.1) then
    write(434,*)b
    allocate(finto(n,n))
    allocate(orto(n,n))
    finto=b
    !    end if
    !  call dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
    call  rgg(nm,n,a,b,alphar,alphai,den,matz,z,info)
    eigen=(alphar)/den
    if (info.ne.0)write(6,*)'info diagol=',info
    write(18,*)eigen
    ! for rgg (if b not positive definite )
    !---------------------------------------------------------------------
    !write(6,*)'ginger' 
    write(28,*)z
    !---------------------------------------------------------------------
    !    if (counter2.eq.1) then

    !-------------------------------------------------------
    !Calculation of the scalar product of the open functions
    !-------------------------------------------------------
    do l=1,n
       do f=1,n 
          orto(l,f)=0.0d0
          do i=1,n
             do j=1,n  
                orto(l,f)=orto(l,f)+z(i,l)*finto(i,j)*z(j,f)
             end do
          end do
          if (orto(l,f).lt. 0.0d0) then
             orto(l,f)=-orto(l,f)
          end if
          !             write(371,*)orto(l,f),l,f
       end do
    end do
    !--------------------------------------------------------
    !Orthogonalization and sorting of the functions
    !--------------------------------------------------------

    !       call ortog(z,finto,n,orto)

    !--------------------------------------------------------
    !Control of the orthogonalization
    !--------------------------------------------------------

    do l=1,n
       do f=1,n 
          orto(l,f)=0.0d0
          do i=1,n
             do j=1,n
                orto(l,f)=orto(l,f)+z(i,l)*finto(i,j)*z(j,f)
             end do
             if (mod(i,2).ne.0) then
             end if
          end do
          if (orto(l,f).lt. 0.0d0) then
             orto(l,f)=-orto(l,f)
          end if
          !             write(373,*)orto(l,f),l,f
       end do
    end do

    if (counter2.eq.1) then

       call sort_shell(eigen,z,finto)


       write(19,*)'eigen',eigen
       write(27,*)z
       allocate(open_vect(n,n))
       open_vect=z
       !------------------------------------------------------------
    else 
       write(23,*)'eigen',eigen,E
       deallocate(orto)
       deallocate(finto)
       allocate(eigen_vector(n,n))
       allocate(eigen_value(n))
       eigen_vector=z
       eigen_value=eigen
       !call cross_section
       deallocate(eigen_vector)
       deallocate(eigen_value)

    end if
    write(21,*)'alphai=',alphai 

    !    a=z

    ! for dsygv (if b positive definite)
    !---------------------------------------------------------------------
    !    write(16,*)'vect',w 
    !    write(15,*)a
    !    open_vect=a 
    !---------------------------------------------------------------------
    deallocate(eigen)
    deallocate(work)
    deallocate(w)   
    deallocate(z)    
    deallocate(vr)
    deallocate(vl)
    deallocate(alphai)
    deallocate(alphar)
    deallocate(den) 
  end subroutine diagol


  subroutine sort_shell(eigen_vect,open_vect,overlap)
    use nrtype
    implicit none 
    real(kind=dbl),dimension(:),INTENT(INOUT) :: eigen_vect
    real(kind=dbl),dimension(:,:),INTENT(INOUT) :: open_vect,overlap 
    real(kind=dbl),dimension(:),allocatable :: prov,provola
    integer(kind=i4b) :: i,j,inc,n,l,max
    real(kind=dbl) :: V
    n=size(eigen_vect)
    allocate(prov(n))
    allocate(provola(n))
    inc=1
    max=n
    do 
       inc=3*inc+1
       if (inc > n) exit
    end do
    do 
       inc=inc/3
       do i=inc+1,n
          v=eigen_vect(i)
          do l=1,max
             prov(l)=open_vect(l,i)
             provola(l)=overlap(l,i)
          end do
          j=i
          do 
             if (eigen_vect(j-inc) <= v) exit
             eigen_vect(j)=eigen_vect(j-inc)
             do l=1,max
                open_vect(l,j)=open_vect(l,j-inc)
                overlap(l,j)=overlap(l,j-inc)
             end do
             j=j-inc
             if(j <= inc) exit
          end do
          eigen_vect(j)=v
          do l=1,max
             open_vect(l,j)=prov(l)
             overlap(l,j)=provola(l)
          end do
       end do
       if (inc <= 1) exit
    end do
    deallocate(prov)
    deallocate(provola)
  end subroutine sort_shell



  subroutine cross_section(J,Jp,N,Np)
    use nrtype, only : dbl,i4b,nep_e,Pi,dpc
    use Brep,only : R0,epsilon_custom
    use Open_information,only : kappa,E,max_open,lmax
    use control , only : cont,molecule,option_Coulomb
    implicit none
    integer(kind=i4b) :: m,p,k,ll,multiplicity,princ_quant_num,channel,xnu,i
    real(kind=dbl) :: x,rj,ry,rjp,ryp,maxo,Bi,Aa,E_ryd,W
    real(kind=dbl),dimension(max_open):: J,Jp,N,Np
    real(kind=dbl),dimension(:),allocatable ::pre_J,pre_Jp,pre_N,pre_Np
    real(kind=dbl),parameter :: ac=1.0d-14
    complex(kind=dpc) :: lambda
    !--------------------------------
    ! Subroutine that calls the ones that calculate the coulomb anr bessel
    ! functions for the matching at the sphere's boundary  
    !--------------------------------
    allocate(pre_J(lmax))
    allocate(pre_Jp(lmax))
    allocate(pre_N(lmax))
    allocate(pre_Np(lmax))
    !write(6,*)'molecule=',molecule
    do m=1,lmax
       xnu = m-1
       if (molecule.eq.'neutral') then
          kappa=sqrt(2.0d0*E)
          x=R0*kappa
          call sphbes(xnu, x,rj,ry,rjp,ryp)
          pre_J(m)=rj
          pre_Jp(m)=kappa*rjp
          pre_N(m)=ry
          pre_Np(m)=kappa*ryp 
       else
          !          write(6,*)'Z=',Z
          !write(6,*)'inside coulomb',Z,R0,E,xnu
          x=R0*Z
          E_ryd=2.0d0*E/Z**2
          kappa=sqrt(E_ryd)
          !--------------------------------------------------------------------------------------------
!!$ Coulomb functions calculation
          ! 1) Option1: Seaton's code 
!          if (option_Coulomb.eq.'Seaton') then
!             call coulomb(cont,xnu,E_ryd,x,ac,rj,ry,rjp,ryp)
             ! 2) Option2: Barnett's code
!          else
             lambda=cmplx(dble(xnu),0.d0)
             rj=0.d0;ry=0.d0;rjp=0.d0;ryp=0.d0
             call call_coulomb_complex(E,lambda,Z,R0,rj,ry,rjp,ryp)
!          end if
          !---------------------------------------------------------------------------------------------
          !write(306,*)E,lambda,Z,R0,(rj),(ry),(rjp),(ryp),lambda
          !*********
          !write(6,*)'Wronskian=',rj*ryp-ry*rjp
          Aa=1.d0
          do i=0,m-1
             Aa=Aa*(1.0d0+E_ryd*1.d0*(i**2))
          end do

          Bi=Aa/(1- exp((-2*Pi/kappa)))
          !          write(6,*)Bi,E
          rj=1.d0/ sqrt(2.d0)*sqrt(Bi)*rj
!*sqrt(2.d0) ! multiply by 2^0.5 to have Wronskian =2/Pi
          rjp=1.d0/ sqrt(2.d0)*sqrt(Bi)*rjp
!*sqrt(2.d0) ! multiply by 2^0.5 to have Wronskian =2/Pi

          ry=1.d0/ sqrt(2.d0*Bi)*ry
!*sqrt(2.d0)   ! multiply by 2^0.5 to have Wronskian =2/Pi
          ryp=1.d0/ sqrt(2.d0*Bi)*ryp
!*sqrt(2.d0)  ! multiply by 2^0.5 to have Wronskian =2/Pi
          if ((abs(rj*ryp-ry*rjp).gt.(1.d0/Pi+epsilon_custom)))write(6,*)'Wronskian error=',rj*ryp-ry*rjp

          !          write(6,*)m-1,'fl=',rj,'gl=',ry,'fjp=',rjp,'fyp=',ryp

          !*COMM*** 
          ! Transformations to obtain Coulomb functions to match the actual wavefunction,
          ! not u=psi*r 

          pre_J(m) = rj/R0
          pre_Jp(m) = Z*rjp/R0-rj/(R0*R0)
          pre_N(m) = -(ry/R0)
          pre_Np(m) = -(Z*ryp/R0-ry/(R0*R0))

          ! *COMM*** 
          ! Sign is changed to the irregular coulomb function because of a different
          !   convention respect to the one used by Seaton (MJ Seaton Rep. prog. phys.
          !   46,167 (1983)) 


          ! Calculate wronskians

          W=rj*ryp-ry*rjp
          !*********
       endif
    end do

    channel=0
    !    princ_quant_num=lmax
    !    channel=0
    !    do k=1,princ_quant_num
    do ll=0,lmax-1
       multiplicity=2*ll+1
       do m=1,multiplicity
          channel=channel+1
          J(channel)=pre_J(ll+1)
          Jp(channel)=pre_Jp(ll+1)
          N(channel)=pre_N(ll+1)
          Np(channel)=pre_Np(ll+1)
       end do
    end do




    deallocate(pre_J)
    deallocate(pre_Jp)
    deallocate(pre_N)
    deallocate(pre_Np)


    end subroutine cross_section


  end module Solve
