!!!*************************************************************
! 文件/File: logdvp.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: logdvp.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

subroutine logdvp (n, n_open, tr, ti, myrank, nproc)
  ! routine to solve the close coupled equations using an
  ! improved log derivative algorithm. the diagonal of the
  ! coupling matrix evaluated at the midpoint of each sector
  ! is used as a reference potential for the sector.  
  use precision
  use shared 
  use molecular_parameters
  use channels, eint=>local_energy 
  use long_range_mod
  use var_ph_module

  implicit none
  integer :: n, n_open
  integer, intent(in)  :: myrank, nproc  ! mpi variables

  real(kind=wp), dimension(n,n) :: y
  
  ! n is the total number of channels
  ! (on return) y is log-derivative matrix Z=h*Y
  ! myrank and nproc are needed to correctly handle the parallel enviroment
  !   in the case a parallel version of main() is built
  
  ! internal variables: 
  ! u        -> potential matrix (output form vpot)
  ! diag     -> reference diagonal potential (output form vpot)
  ! eshift   -> (ecoll - k^2)*2*mu

  real(kind=wp), dimension(n,n) :: u
  real(kind=wp), dimension(n) :: y14, y23, eshift, diag
  real(kind=wp) :: esave, h, d1, d2, d4, half, sign, hi
  real(kind=wp) :: wref, r, arg, tn, th
  integer :: istart, nodes, nsave, ir, i, j, kstep, kount, unit  
  real(kind=wp), intent(out), dimension(n_open,n_open) :: tr, ti

  ! Local variables
  real(kind=wp), allocatable, dimension(:,:) :: k, k_open
  real(kind=wp) ::  r_airy, hmin, rout, delta
  integer       :: ii
  ! statistic of var_ph
  integer       :: NIV, NLDR, NCH, nok, nbad

  allocate (k(n,n))
  allocate(k_open(n_open,n_open))
  beszero=.true.

  !unit on which we write debug information
  unit=myrank+20
  if (nproc==1) unit=6  ! std output

999 format(t30,a)
  if (debug>=1) then
     write(unit,999)'+-----------------+'
     write(unit,999)'| entered Log-der |'
     write(unit,999)'+-----------------+'
  endif
  istart=0
  nodes=0
  esave=reduced_e
  ! here we put k^2 diagonal matrix into the equations.
  do i=1,n
     eshift(i)=eint(i)-reduced_e
  enddo

  ! this version uses a constant step size through the
  ! integration range, with nsteps steps between rbegin and rswitch.
  h=half_step_size
  d1=h*h/3.0_wp
  d2=2.0_wp*d1
  d4=-d1/16.0_wp
  half=0.5_wp*h
  nsave=0
  r=rbegin
  ir=1
  ! r  -> is the value of r
  ! n  -> is the dimension of the scattering problem
  ! on output u    -> is the matrix for propagator
  ! on output diag -> is the diagonal part of potential
  call vpot(r, n, u, diag)
  u=d1*u
  ! istart=0/1 means that initial log derivative matrix isn't/is
  ! already in y.  default is 0.
  if(istart/=1) then
     sign=1.0_wp
     if(rend<rbegin) sign=-1.0_wp
     y=0.0_wp    ! initialize y
     do j=1,n
        wref=diag(j)+eshift(j)
        y(j,j)=sign*1.e+30_wp
        if(wref>0.0_wp) y(j,j)=sign*sqrt(wref)
     enddo
  endif
  y=h*y+u 
  ! main agorithm is unchanged
  cycle_over_r: do kstep=1,nsteps
     r=r+h
     ir=ir+1
     call vpot(r, n, u, diag)
     do j=1,n
        do i=j,n
           u(i,j)=d4*u(i,j)
        enddo
     enddo
     do i=1,n
        u(i,i)=0.125_wp
     enddo
     call syminv(u,n,n)
     do i=1,n
        u(i,i)=u(i,i)-8.0_wp
     enddo
     do i=1,n
        wref=diag(i)+eshift(i)
        arg=half*sqrt(abs(wref))
        if (wref<0.0_wp) then
           tn=tan(arg)
           y14(i)=arg/tn-arg*tn
           y23(i)=arg/tn+arg*tn
        else
           th=tanh(arg)
           y14(i)=arg/th+arg*th
           y23(i)=arg/th-arg*th
        endif
        u(i,i)=u(i,i)+2.0_wp*y14(i)
        y14(i)=y14(i)-d1*diag(i)
        y14(i)=max(y14(i),0.0_wp)
        y(i,i)=y(i,i)+y14(i)
     enddo
     call syminv(y,n,n)
     do j=1,n
        do i=j,n
           y(i,j)=u(i,j)-y23(i)*y(i,j)*y23(j)
        enddo
     enddo
     call syminv(y,n,n)
     r=r+h
     ir=ir+1
     call vpot(r, n, u, diag)
     if (kstep.eq.nsteps) d2=d1
     do j=1,n
        do i=j,n
           u(i,j)=d2*u(i,j)
        enddo
     enddo
     do j=1,n
        do i=j,n
           y(i,j)=u(i,j)-y23(i)*y(i,j)*y23(j)
        enddo
     enddo
     do i=1,n
        y(i,i)=y(i,i)+y14(i)
     enddo
  enddo cycle_over_r
  hi=1.d0/h  ! next cycle forms the Y-matrix=Z/h
  do j=1,n
     do i=j,n
        y(i,j)=hi*y(i,j)
        y(j,i)=y(i,j)
     enddo
  enddo
  if (debug>=1) then
     write(unit,999)'Exiting Log-der'
  endif
  istart=0  

  ! This part handles the long_range part of the propagation
  ! depending on the value of select_propagator variable.
  ! If select_propagator='PM' then pure Manolopulos has been chosen
  !   and in this case only log-derivative with constant stepsize is used.
  !   The log_der() subroutine is called and used to propagate solution up
  !   to rswitch (rend is ignored in this case). Here the Y matrix is 
  !   transformed in K matrix using ytok() and then in T-matrix 
  !   using ktos().
  !   Note that in this and in the next case the K-matrix is derived 
  ! only for the Open-Open block. 
  ! Manolopulos and Variable Phase scheme has been 
  !   chosen and we get the K matrix at rswitch using ytok_aug (augmented)
  !   which retursns the "complete" K-matrix containing also Close-Open 
  !   and Closed-Closed blocks. We then call the Variable Phase propagator
  !   up to rend and transform then K-matrix using ktos().

     ! +----------------------------------------------+ 
     ! | Manolopulos+Variable phase: Y->K  vp()  K->T |
     ! +----------------------------------------------+   
     beszero=.true.

     if (debug>=1) then
        write(unit,999)'Calculating K-matrix'
     endif
     ! calculating complete K-matrix

     wave_vector=sqrt(abs(reduced_e-eint))

     call ytokaug(wave_vector, l_local, n, n_open, y, k, rswitch, beszero)

     if (debug>=1) then
        write(unit,999)'Calling variable phase propagator'
     endif
     ! call variable phase propagator

     call varp_rk(k, k_open ,n ,rswitch, rout,  &
          nok, nbad, n_open, NIV, NLDR, NCH, rend, myrank, nproc) 
          !these last are returned values
 
     if (debug>=1) then
        write(unit,*) ' '
        write(unit,1001) 'end of variable phase propagation at rend= ', rout
        if(onlymod) then
           write(unit,1000) 'user controlled end of the propagation'
        else
           write(unit,1000) 'automatic end of the propagation'
        endif
        write(unit,*) ' '
        write(unit,1000) 'statistics of variable phase-propagation'
        write(unit,1000) '========================================'
        write(unit,1002) 'nok=',   nok
        write(unit,1002) 'nbad=',  nbad 
        write(unit,1002) 'number of inversions=',niv
        write(unit,1002) 'number of calls to log-der=',nldr
        write(unit,1002) 'final number of channels retained=',nch
        write(unit,1002) 'number of open channels=',n_open
        write(unit,*) ' '                          
     endif

     call ktos(k_open,tr,ti,n_open)

     if (debug>=2)then
        write(unit,*)' '
        write(unit,999)'+-----------+'
        write(unit,999)'| s-matrix: |'
        write(unit,999)'+-----------+'
        write(unit,*)' '
        write(unit,1051)' i',' ii',' re(s)',' imm(s)'
     endif
     do i=1,n_open
        do j=1,n_open
           if (debug>=2)then
              if (i.le.j)then
                 write(unit,1052)i,j,tr(i,j),ti(i,j) 
              endif
           endif
           ! now tranform S into T
           delta=0.0_wp
           if (i.eq.j) delta=1.0d0
           tr(i,j)=tr(i,j)-delta
           ti(i,j)=ti(i,j)
        enddo
     enddo
     deallocate(k)
     deallocate(k_open)
     return
!999 format(t30,a)
1000 format(t35,a)
1001 format(t35,a,f15.3)
1002 format(t35,a,i6)
1051 format(t30,a,t40,a,t50,a,t65,a)
1052 format(t30,i3,t40,i3,t50,e12.6e2,t65,e12.6e2)
end subroutine logdvp











