!!!*************************************************************
! 文件/File: derivs1.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: derivs1.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

subroutine kder_check(r,y,dydx,nvar,n,nopen,ninit,nch_on_shell,flagch)
  use precision
  use channels, l=>l_local, wvec=>wave_vector
  use var_ph_module
  implicit none
  !--------------------------------------------------------------
  !     this routine computes the derivative of the 
  !     scaled-augmented k-matrix to be used in 
  !     the variable phase propagator
  !     
  !     on input:
  !     --------
  !     
  !     - r is the independent variable (scattering coordinate)
  !     - y is an array of dimension m=n*(n+1)/2
  !       (size of the differential problem)
  !       that contains the upper triangle of the k
  !       matrix (k(m,m))
  !       y(1)  y(2)  y(4)  y(7)   ...
  !             y(3)  y(5)  y(8)   ...
  !                   y(6)  y(9)   ...
  !                         y(10)  ...
  !     
  !       (this choiche allows in a simple way the 
  !       reduction of dimensionality)
  !     - nvar is the dimensionality of the ode problem
  !     - n is the actual number of channels
  !     - nopen is the number of open channels
  !     - ninit is the original number of channels
  !       (the first n of which are now used)
  !     
  !     on output:
  !     ---------- 
  !     
  !     - dydx contains the derivatives of the k-matrix elements
  !     
  !------------------------------------------------------------------
  !     input/output arrays
  integer, intent(in) :: nvar, nopen, ninit, n, nch_on_shell
  real(kind=wp), intent(inout), dimension(n) :: y, dydx
  real(kind=wp), intent(in) :: r
  logical, intent(out) :: flagch

  
  !     local variables
  real(kind=wp), dimension(n,n) :: yk, u, u1, a, aclo, psi, b
  real(kind=wp), dimension(n)   :: sj, sn, diag
  real(kind=wp) :: dw, darg, rootdw, rj, rn, factor, rjp, rnp, rsq
  integer :: nm, i, ii, j, nclose
  
  !------------------------------------------------------------------
  ii=1
  do j=1,n
     do i=1,(j-1)
        yk(i,j)=y(ii)
        yk(j,i)=y(ii)
        ii=ii+1
     enddo
     yk(j,j)=y(ii)
     ii=ii+1
  enddo
  !     
  !     compute arguments of bessel functions
  !     
  nclose = n - nopen
  !
  if(beszero) then
     do i = 1,nopen
        dw = wvec(i)
        darg = dw*r
        rootdw = sqrt(dw)
        rj=sin(darg)
        rn=-cos(darg)
        if(realk) then
           sj(i) = rj/rootdw
           sn(i) = rn/rootdw
        else
           sj(i) = rn/rootdw
           sn(i) = rj/rootdw
        endif
     enddo
     if(nclose.gt.0) then
        factor=sqrt(2.0_wp/acos(-1.0_wp))
        do i = (nopen+1),n
           dw = wvec(i)
           darg = dw*r
           rj=factor*(1.0_wp-exp(-2.0_wp*darg))/2.0_wp
           rn=-1.0_wp/factor
           rootdw = sqrt(dw)
           if(realk) then
              sj(i) = rj/rootdw
              sn(i) = rn/rootdw
           else
              sj(i) = rn/rootdw
              sn(i) = rj/rootdw
           endif
        enddo
     endif
     !
  else
     !
     do i = 1,nopen
        dw = wvec(i)
        darg = dw*r
        call besopen(darg, l(i), rj, rn, rjp, rnp)
        rootdw = sqrt(dw)
        if(realk) then
           sj(i) = rj/rootdw
           sn(i) = rn/rootdw
        else
           sj(i) = rn/rootdw
           sn(i) = rj/rootdw
        endif
     enddo
     if(nclose.gt.0) then
        do i = (nopen+1),n
           dw = wvec(i)
           darg = dw*r
           call besclosed(darg, l(i), rj, rn, rjp, rnp)
           rootdw = sqrt(dw)
           if(realk) then
              sj(i) = rj/rootdw
              sn(i) = rn/rootdw
           else
              sj(i) = rn/rootdw
              sn(i) = rj/rootdw
           endif
        enddo
     endif
  endif
  !     
  !     compute psi=j-n*k 
  !     
  do j=1,n
     do i=1,n
        psi(i,j)=-yk(i,j)*sn(i)
        u1(i,j)=0.0_wp
     enddo
     psi(j,j)=psi(j,j)+sj(j)
  enddo
  
  !     
  !     we call vpot to get the potential
  !     matrix with the actual size of the problem
  !           
  call vpot(r, n, u, diag)

  ! vpotlr returns with cetrifugal barriers, 
  ! subtract them if non-modified propagation is requested
  if(.not.beszero) then
     rsq=1.0_wp/r/r
     do i=1,n
        u(i,i)=u(i,i)-real(l(i)*(l(i)+1),kind=wp)*rsq
     enddo
  endif
  
  !     
  !     compute psi_transpose*v*psi
  !     

  flagch=.false.
  close: if(realk.and.nclose.gt.0) then
     !
     ! this is for true matrix propagation
     ! in this case we can remove weakly coupled closed channels
     ! nch_on_shell is the number of "local" channels
     ! in corresponding to the highest asymptotic channel
     ! see also calling routine var_ph
     !
     flagch=.true.
     do j=1,nch_on_shell
        do i=1,n
           u1(i,n+1-j)=u(i,n+1-j)
           u1(n+1-j,i)=u1(i,n+1-j)
           u(i,n+1-j)=0.0_wp
           u(n+1-j,i)=0.0_wp
        enddo
     enddo
     
     call dgemm('n','n',n,n,n,1.0_wp,u,  n,psi,n, 0.0_wp, b,n)
     call dgemm('t','n',n,n,n,1.0_wp,psi,n,b,  n, 0.0_wp, a,n)
     call dgemm('n','n',n,n,n,1.0_wp,u1,  n,psi,n, 0.0_wp, b,n)
     call dgemm('t','n',n,n,n,1.0_wp,psi,n,b,  n, 0.0_wp, aclo,n)
     
     do j=1,n
        do i=1,j
           if(j.le.nopen) then
              if(flagch) then 

                 ! if(abs(aclo(i,j)/a(i,j)).lt.(sigch*nch_on_shell)) then

                 ! this is a small patch from Bodo: the line above gives 
                 ! floating invalid whenever a(i,j)=0. When this happens, it
                 ! seems that also aclo(i,j)=0. So set flagch when they are both zero.

                 if(abs(aclo(i,j)).le.(abs(a(i,j))*sigch*nch_on_shell)) then
                    flagch=.true.
                 else
                    flagch=.false. 
                 endif
              endif
           endif
           a(i,j)=a(i,j)+aclo(i,j)  
        enddo
     enddo
 900  format(2i4,3(f18.7,x))
     !     
     !     this is the term that arises from the use 
     !     of scaled modified bessel functions.
     !     note that it contribuites only to the open-closed
     !     and closed-closed block of the k' matrix
     !     
     do j=(nopen+1),n
        do i=1,j
           a(i,j)=a(i,j)+yk(i,j)*wvec(j)
        enddo
        do i=(nopen+1),j
           a(i,j)=a(i,j)+yk(i,j)*wvec(i)
        enddo
     enddo
  else
     !
     ! this is for inverse k matrix propagation
     ! (it is not really used..see var_ph)
     !
     call dgemm('n','n',n,n,n,1.0_wp,u,  n,psi,n, 0.0_wp, b,n)
     call dgemm('t','n',n,n,n,1.0_wp,psi,n,b,  n, 0.0_wp, a,n)
     !     
     !     this is the term that arises from the use 
     !     of scaled modified bessel functions.
     !     note that it contribuites only to the open-closed
     !     and closed-closed block of the k' matrix
     !     
     do j=(nopen+1),n
        do i=1,j
           a(i,j)=a(i,j)+yk(i,j)*wvec(j)
        enddo
        do i=(nopen+1),j
           a(i,j)=a(i,j)+yk(i,j)*wvec(i)
        enddo
     enddo
  endif close
  !     
  !     fill dydx=-(psi_transpose*v*psi+scale_factor)
  !     
  ii=1
  if(realk) then
     do j=1,n
        do i=1,(j-1)
           dydx(ii)=-a(i,j)
           ii=ii+1
        enddo
        dydx(ii)=-a(j,j)
        ii=ii+1
     enddo
  else
     do j=1,n
        do i=1,(j-1)
           dydx(ii)=a(i,j)
           ii=ii+1
        enddo
        dydx(ii)=a(j,j)
        ii=ii+1
     enddo
  endif
  return  
end subroutine kder_check













