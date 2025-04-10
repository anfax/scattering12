!!!*************************************************************
! 文件/File: daprop_lr.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: daprop_lr.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

subroutine daprop_lr (n, ninit, n_open, yk, rx) 
  ! routine to solve the close coupled equations using an
  ! improved log derivative algorithm. the diagonal of the
  ! coupling matrix evaluated at the midpoint of each sector
  ! is used as a reference potential for the sector.  
  use precision
  use shared 
  use molecular_parameters
  use channels, l=>l_local, eint=>local_energy
  use var_ph_module

  implicit none
  integer, intent(IN) :: n, n_open, ninit
  real(kind=wp), dimension(ninit,ninit), intent(inout) :: yk
  real(kind=wp), intent(inout) :: rx

  ! n is the total number of channels
  ! n_open are the open channels
  ! n_lambda is the number of symmetries used
  ! eint() is e_channel*cz=e_channel*2*mu/hbar^2
  ! v is the potential_coefficient(lambda,n,n)

  ! internal variables: 
  ! u        -> potential matrix (output form vpot)
  ! y        -> log-der matrix
  ! diag     -> reference diagonal potential (output form vpot)
  ! eshift   -> (ecoll - k^2)*2*mu

  real(kind=wp), dimension(n,n) :: u, y, yk1
  real(kind=wp), dimension(n) :: y14, y23, eshift, diag
  integer :: istart, nodes, nsave, ir, i, j, kstep_lr, kount, unit, ii, nsteps_lr
  real(kind=wp) :: esave, h, d1, d2, d4, half, sign
  real(kind=wp) :: wref, arg, tn, th, hi, delta
  real(kind=wp) :: hmin, rout, sum
  
  nodes=0
  esave=reduced_e
  do i=1,n
     eshift(i)=eint(i)-reduced_e
     do j=1,n
        yk1(j,i)=yk(j,i)
     enddo
  enddo
 
  nsteps_lr=10
  h=step_size
  d1=h*h/3.0_wp
  d2=2.0_wp*d1
  d4=-d1/16.0_wp
  half=0.5_wp*h
  nsave=0
  ir=1

  call vpot(rx, n, u, diag)

  u=d1*u
 
  call KTOY(wave_vector, l, n, n_open, yk1, y, rx, beszero) 

  y=h*y+u 
  cycle_over_r: do kstep_lr=1,nsteps_lr
     rx=rx+h
     ir=ir+1
     call vpot(rx, n, u, diag)
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
     rx=rx+h
     ir=ir+1

     call vpot(rx, n, u, diag)
     if (kstep_lr.eq.nsteps_lr) d2=d1
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

  hi=1.d0/h
  do j=1,n
     do i=j,n
        y(i,j)=hi*y(i,j)
        y(j,i)=y(i,j)
     enddo
  enddo

  call YTOKAUG(wave_vector, l, n, n_open, y, yk1, rx, beszero)
  
  do j=1,n
     do i=1,n
        yk(i,j)=yk1(i,j)
     enddo
  enddo 

end subroutine daprop_lr

      

