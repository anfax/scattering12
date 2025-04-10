!!!*************************************************************
! 文件/File: var_ph.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: var_ph.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

subroutine varp_rk(yk, ykout, ninit, rbegin, rout, &
     nok, nbad, nopen, niv, nldr, nch, rvarp, myrank, nproc)
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !     
  !     modified variable phase propagator for long-range potentials.
  !     it solves the scaled and modified  k-matrix variable 
  !     phase equation using a cash-karp embedded runge-kutta algorithm
  !     (see chapter 16 of "numerical recipes in fortran 77: 
  !     the art of scientific computing" copyright(c) 1986-1992 by 
  !     cambridge university press, also at http://www.nr.com)
  !
  !     if the required stepsize is lower than a chosen value (hmin) the
  !     method automatically propagates the inverse of the k matrix
  !     solving the corresponding vp equations. a safety log-derivative 
  !     call is inserted to avoid infinite loops near a singular point.
  !
  !     the variable phase equations are scaled in the sense that
  !     the closed channel basis functions are exponentially damped
  !     to avoid numerical overflows/underflows; this results in slightly
  !     modified equations, accordingly to 
  !     "the variable phase method in multichannel electron atom scattering"
  !     m. le dourneuf and vo ky lan.
  !
  !     more important, the equations are modified in the sense their
  !     computational cost is lowered by avoiding the calculation of 
  !     the bessel functions for open channels for the main part of 
  !     propagation. this has been reached by using an effective 
  !     potential that includes the centrifugal potential in the open 
  !     channel space.
  !     moreover the number of channels effectively used in the propagation
  !     is reduced "on the fly" whenever is possible. it is accomplished 
  !     by computing the contribution of the highest energy closed channel
  !     to the open-open block of the k-matrix and by comparing it with
  !     the other contributions.
  !     
  !     
  !
  !     on input:
  !     --------
  ! 
  !     - yk is the starting k-matrix. it comes from the log-derivative
  !       matrix of the inner problem through a suitable ytok transformation.
  !       it fills the y arrays that contains the symmetry unique matrix
  !       elements. the filling is as follows:
  !     
  !                y(1) y(2) y(4) ..
  !                     y(3) y(5) ..
  !                          y(6) ..
  !                               ..
  !       this choice allows a simple implementation 
  !       of the channel reduction procedure
  !     
  !     - ninit is the dimension of the initial k-matrix (ninit*ninit)
  !     
  !     - rbegin is the starting point of integration
  !     
  !     * eps (by module var_ph_module) is a tollerance 
  !       variable that sets the required accuracy of the i-th equation as
  !     
  !               delta=eps*(abs(y(i))+abs(h*dydx(i)))
  !     
  !     * h is the starting stepsize (from module shared)
  !     
  !     * hmin is the minimum allowed stepsize (can be zero).  
  !       it should be set equal to the log-der step. 
  !       when the adaptive step procedure leads to a value lower
  !       than this one the propagator switch to the k^-1 propagation.
  !       (it may happen that the optimum stepsize for the k^-1 
  !       propagation is also lower than hmin; in this case some
  !       logderivative steps are used to overcome the singular 
  !       region) 
  !     
  !     * soglia (from module var_ph_module) is the treshold value 
  !       below which the current step 
  !       is a good candidate for ending the propagation
  !
  !                    abs(dydx(i)).lt.soglia          for all i 
  !
  !       is the relation checked.
  !       (nendmax fixes the number of consecutive steps
  !       needed to end the propagation)
  !     
  !     * nopen is the number of open channels
  !     
  !     * beszero (from module var_ph_module) is a logical variable 
  !       that defines the modified 
  !       algorithm. on input is set to .true. by the calling routine
  !       which means that the modified algorithm (i.e. no bessel
  !       computation for open channels) is being used. it is set
  !       to .false. when r.gt.varp (corresondingly a ktok transformation 
  !       is performed)
  !
  !     * sigch (from module var_ph_module) is the threshold for the channel 
  !       reduciton procedure (see derivs1)
  !
  !     * onlymod(from module var_ph_module)   
  !                 .true.  means that the propagation is only 
  !                         of modified_k kind
  !                 .false. means that automatic ending procedure works
  !                         (it can only work with ordinary vp propagation)
  !
  !     on output:
  !     --------- 
  !       
  !     - ykout is the (physical) k-matrix, i.e. the open-open block of the 
  !       k-matrix at the end of the propagation
  !     
  !     - rout is the ending point of integration (an authomatic procedure 
  !       for such value has been implemented)
  !
  !     - nok is the number of good ode steps
  !
  !     - nbad is number of failed (and retried) steps of ode integration
  !
  !     - niv is the number of matrix inversions (ktok-1 and viceversa)
  !
  !     - nldr is the number of "safety" log-derivative steps
  !
  !     - nch is the final number of channels
  !
  !    
  !
  !     calls to:
  !     --------
  !
  !     -kder_check:   the routine that computes the derivative 
  !                 of the k-matrix (or the derivative
  !                 of the inverse k-matrix if it is the case). 
  !                 it differs sliglthy from the routine derivs
  !                 called by the ode integrator because it 
  !                 computes the contributions of the last 
  !                 closed channel (we assume that they are energy-sorted)
  !                 to the open-open block (of k or k-1) and compares this
  !                 values with the remaining contributions:
  ! 
  !                    dk/dx(i,j)= dk/dx(i,j)_from_last_channel+
  !                                dk/dx(i,j)_from_other_channels
  !                 (i,j: open)
  !                 if (dkdx_last_channel).lt.(sigch*dkdx_open_channel)
  !                 then the last channel is a good candidate
  !                 for the channel reduction procedure.
  !
  !     -rkqs:      the routine that makes a "quality controlled"
  !                 runge-kutta step using the cash-karp prescriptions
  !                 (numerical recipes).it calls
  !                 - rkck which actually performs the integration, 
  !                      by calling
  !                     -derivs which computes the "ordinary" 
  !                      derivative of k 
  !                      or k-1 by calling 
  !                        - vpot_lr (potential evaluation)
  !                        - bessik, bessij (bessel evaluation) 
  !
  !     -daprop_lr: an adapted constant reference log-derivative propagator
  !                 that calls vpot_lr for potential evaluation. it performs
  !                 also the initial ktoy step and the final ytok step.
  !
  !     -syminv:    to do  symmetric matrix inversions
  !     
  !     -ktok:     perform a modified_k-to-k matrix transformation.
  !                 for simplicity it performs a modified_k-to-y step
  !                 and a y-to-k step. it is needed at the 
  !                 end of the modified propagation. 
  !
  !     internal variables:
  !     -----------------
  !      
  !     - realk     : .true. means that the k matrix is being propagated.
  !                   .false. means that the inverse k matrix is being 
  !                           propagated.
  !                    default is .true.
  !                    compatible with beszero=t/f flag.
  !
  !     - flach     :  output from derivs1 (default is .false.) .
  !                    .true. means that the last channel is    
  !                           a good candidate
  !                           to be eliminated. if it remains true for 5
  !                           consecutive steps the channel is eliminated
  !                    compatible with beszero=t/f flag.
  !
  !
  !     - ftresh    : .true. means that all the open-open k-matrix elements
  !                          satisfy the ending condition. 
  !                          if it remains true
  !                          for nendmax consecutive steps 
  !                          the propagation is ended.
  !                    it acts only when beszero=.false. 
  !                    default is .false.
  !
  !     - safety    : integration check. it is set to .false. when the step 
  !                   becomes lower than hmin and a k-to-invk transformation
  !                   (or viceversa) is needed. it is set to .true. if the 
  !                   successive step is successful, otherwise it enables
  !                   the call to the log-der propagator. default is .true. 
  !
  !     - first     : .true. is the default. It is set to false when
  !                   when flagch is set to .true. for the first time 
  !                   at a certain xsav (i.e. locally). 
  !                   It allows to measure the interval dxsav
  !                   from xsav in which flagch remains true.
  !                   
  !     - nvar      : number of od equations
  !
  !     - nvaropen  : number of od equations for open-open matrix elements
  !
  !     - nend      : counter for the ending procedure (see ftrash)
  !
  !     - nch_on_shell: the number of the (closed) channels in the 
  !                     outest shell.
  !                     Here "shell" means the ensemble of channels 
  !                     with the same wavevector. 
  !  
  !     - xsav      : it the starting point from the checking of the channel
  !                   reduction procedure of the last nch_on_shell channels.
  !                   It is the x value for which the requirments on the 
  !                   outest shell are satisfied.
  !
  !     - dxsav     : it is the length of the interval for 
  !                   which the requirments
  !                   on the outest shell are satisfied. To reduce 
  !                   the channels
  !                   it must be longer than a prescribed value 
  !                   (delr_of_shell)
  !
  !     - delr_of_shell: the minimum interval length for which 
  !                   flagch must remains
  !                   true to reduce the numebr of channels. 
  !                   It is an energy dependent
  !                   quantity set to const/(wavevector of the shell)
  !
  !      
  !---------------------------------------------------------------------
  use precision
  use shared, rdummy=>rbegin
  use var_ph_module
  use channels, l=>l_local, nvnumbers=>n_local, eint=>local_energy, &
       wvec=>wave_vector
  
  implicit none
  integer, intent(in)  :: ninit, nopen, myrank, nproc
  integer, intent(out) :: nok, nbad, niv, nldr, nch
  real(kind=wp), dimension(ninit,ninit), intent(INOUT) :: yk
  real(kind=wp), dimension(nopen,nopen), intent(OUT) :: ykout
  real(kind=wp), intent(in)  :: rbegin, rvarp
  real(kind=wp), intent(out) :: rout
  
  ! local variables
  logical flagch, ftresh, safety, first
  integer nvaropen, nvar, i, ii, knot, nend, j, unit, nstp, nch_on_shell
  real(kind=wp) :: tresh, h, hdid, hnext, x
  real(kind=wp) :: hmin, dxsav, xsav, delr_of_shell
  real(kind=wp), dimension(ninit,ninit):: yk1
  real(kind=wp), dimension(ninit*(ninit+1)/2) :: dydx, y, yscal

  ! some precision parameters
  real(kind=wp), parameter :: tiny=1.e-30
  integer, parameter :: nendmax=5
  unit=myrank+20
  if (nproc==1) unit=6  ! std output

  ! set default values
  h=half_step_size
  hmin=half_step_size
  realk=.true.
  safety=.true.
  ftresh=.false.
  nvar=ninit*(ninit+1)/2
  nch=ninit
  nend=0
  niv=0
  nldr=0
  nok=0 
  nbad=0 
  !determine the number of channels in the outest energy level
  !and set the checking interval for the reduction
  nch_on_shell=0
  if((nch-nopen).gt.0)then
     nch_on_shell=1
     i=1
     do
        if(nch>i.and.wvec(nch-i)==wvec(nch)) then
           nch_on_shell=nch_on_shell+1
           i=i+1
        else
           exit
        endif
     enddo
  endif
  delr_of_shell=0.1_wp/wvec(nch)
  
  nvaropen=nopen*(nopen+1)/2
  x=rbegin 

999 format(t35,a)
1999 format(a,i3,t10,a,f10.3)
1000 format(t35,a,f10.3)
1002 format(t35,a,E10.5)
1001 format(t35,a,i5)
1991 format(a,i3,t10,a,i6)  
  
  if (debug==0) then
     if(beszero)then
        write(unit,1999)'RK-',myrank,'-> Called modified variable phase propagator at r=', & 
             x*au_to_angst
     else
        write(unit,1999)'RK-',myrank,'-> Called variable phase propagator at r=', & 
             x*au_to_angst
     endif
  endif
  if (debug>=1) then
     write(unit,999)'+---------------------------+'
     write(unit,999)'| Variable-phase propagator |'
     write(unit,999)'+---------------------------+'
     if(beszero)then
        write(unit,1000)'called modified variable phase propagator at r=', & 
             x*au_to_angst
     else
        write(unit,1000)'called variable phase propagator at r=', & 
             x*au_to_angst
     endif
     write(unit,1001)'with number of channels               ', ninit
     write(unit,1001)'number of open channels               ',nopen
     write(unit,1000)'minimum allowed step size             ',hmin
     write(unit,1002)'tollerance in od integration          ',eps
     write(unit,1002)'tollerance in channels reduction      ',sigch
     if(.not.onlymod) then 
        write(*,1002)'tollerance in ending condition        ',soglia
     endif
  endif

  ! load the y array with the symmetry unique k-matrix elements
  ii=0
  do j=1,ninit
     do i=1,j
        ii=ii+1
        y(ii)=yk(i,j)
     enddo
  enddo    

  ! start adaptive stepsize integration
  nstp=0
  cycle_over_r: do 
     nstp=nstp+1
     ! end of the modified propagation - 001
     beszero_check: if(beszero) then
        termination_check: if(x.gt.rvarp) then
           ii=0
           do j=1,nch
              do i=1,j
                 ii=ii+1
                 yk(i,j)=y(ii)
                 yk(j,i)=y(ii)
              enddo
           enddo
           if(.not.realk) then
              call syminv(yk,ninit,nch)
              niv=niv+1
              do j=1,nch
                 do i=j,nch
                    yk(j,i)=yk(i,j)
                 enddo
              enddo
              realk=.true.
              safety=.true.
           endif
           call ktok(wvec,l,ninit,nch,nopen,yk,x)
           if(onlymod)then
              do j=1,nopen
                 do i=1,j
                    ykout(i,j)=yk(i,j)
                    ykout(j,i)=ykout(i,j)
                 enddo
              enddo
              rout=x
              ! finishing only modified propagation
              if (debug>=1) then
                 write(unit,1000)'   Modified Var-phase propagation finished at r= ',& 
                      rout*au_to_angst
                 write(unit,1001)'   Number of steps performed: ', nstp
              elseif(debug==0) then
                 write(unit,1999)'RK-',myrank,'   Modified Var-phase propagation finished at r= ', & 
                      rout*au_to_angst
                 write(unit,1991)'RK-',myrank,'   Number of steps performed: ', nstp
              endif
              return
           else
              ii=0
              do j=1,nch
                 do i=1,j
                    ii=ii+1
                    y(ii)=yk(i,j)
                 enddo
              enddo
              beszero=.false.
              h=hmin
           endif
        endif termination_check
     endif beszero_check

     !  "outside" computation of the derivatives - 002
     !  for true k matrix propagation 
     !  it allows also the channel reduction procedure
!    Bodo   
     if(realk)then
        call kder_check(x,y,dydx,nvar,nch,nopen,ninit,nch_on_shell,flagch)
        if(flagch) then
           if(first) then
              dxsav=0.0_wp
              xsav=x
              first=.false.
           else
              dxsav=dxsav+(x-xsav)
              xsav=x
           endif
        else
           first=.true.
        endif
        if (dxsav.gt.delr_of_shell) then
           nch=nch-nch_on_shell
           nvar=nch*(nch+1)/2
           if(debug>=1)then
              write(unit,fmt='(t35,a,f10.2,2x,i5,a,i5)')'var_ph:  channels reduction at r',&
                   x,nch+nch_on_shell,'-->',nch
           endif
           first=.true.
           dxsav=0.0_wp
           delr_of_shell=0.1_wp/wvec(nch)
           nch_on_shell=1
           i=1
           do
              if(nch>i.and.wvec(nch-i)==wvec(nch)) then
                 nch_on_shell=nch_on_shell+1
                 i=i+1
              else
                 exit
              endif
           enddo
        endif
     else
        ! for inverse k matrix propagation we use the faster kder routine
        call kder(x,y,dydx,nvar,nch,nopen,ninit)
     endif
     do  i=1,nvar
        yscal(i)=abs(y(i))+abs(h*dydx(i))+tiny 
     enddo

     
     ! check for automatic end of the propagation - 003 
     if(.not.beszero) then
        ftresh=.true.
        do i=1,nvaropen
           if(ftresh) then
              if(abs(dydx(i)).gt.soglia) ftresh=.false.
           endif
        enddo
        if(ftresh) then
           nend=nend+1
        else
           nend=0
        endif
     endif
     

     ! automatic end of the propagation - 004
     if(nend.eq.nendmax) then
        ii=0
        do j=1,nopen
           do i=1,j
              ii=ii+1
              ykout(i,j)=y(ii)
              ykout(j,i)=y(ii)
           enddo
        enddo
        if(.not.realk) then
           call syminv(ykout,nopen,nopen)
           niv=niv+1
           do j=1,nopen
              do i=j,nopen
                 ykout(j,i)=ykout(i,j)
              enddo
           enddo
        endif
        rout=x
        ! finishing only modified propagation
        if (debug>=1) then
           write(unit,1000)'Var-phase propagation finished at r= ',& 
                rout*au_to_angst
           write(unit,1001)'Number of steps performed: ', nstp
        elseif(debug==0) then
           write(unit,1999)'Var-phase propagation finished at r= ', & 
                rout*au_to_angst
           write(unit,1991)'Number of steps performed: ', nstp
        endif
        return
     endif
     
     
     ! quality controlled runge-kutta step - 006
     call rkqs(y,dydx,nvar,x,h,yscal,hdid,hnext,nch,nopen,ninit)
     if(hdid.eq.h)then
        nok=nok+1
     else
        nbad=nbad+1
     endif

     
     ! handling of possible local resonances - 007
     if(abs(hnext).lt.hmin) then
        ! if hnext.lt.hmin inverse propagation is tried when safety=.true.
        ! otherwise a safety log-d,nch,nopen,ninit)erivative 
        ! propagation is done.
        ! hmin is chosen differently 
        ! according if we are propagating the kmatrix (hmin=half_step_size)
        ! or its inverse (10*hmin): in this way our algorithm 
        ! is unbalanced in favor of the k-matrix propagation
        ! (to perform the channel reduction procedure)
        if(safety) then
           safety=.false. 
           ii=0
           do j=1,nch
              do i=1,j
                 ii=ii+1
                 yk(i,j)=y(ii)
                 yk(j,i)=y(ii)
              enddo
           enddo
           call syminv(yk,ninit,nch)
           niv=niv+1
           do j=1,nch 
              do i=j,nch
                 yk(j,i)=yk(i,j)
              enddo
           enddo
           if(realk) then
              hmin=hmin*10.0_wp
              realk=.false.
           else
              hmin=hmin*0.1_wp
              realk=.true.
           endif
           ii=0
           do j=1,nch
              do i=1,j
                 ii=ii+1
                 y(ii)=yk(i,j)
              enddo
           enddo
        else
           ii=0
           do j=1,nch
              do i=1,j
                 ii=ii+1
                 yk(i,j)=y(ii)
                 yk(j,i)=y(ii)
              enddo
           enddo
           if(.not.realk) then
              call syminv(yk,ninit,nch)
              niv=niv+1
              do j=1,nch
                 do i=j,nch
                    yk(j,i)=yk(i,j)
                 enddo
              enddo
              hmin=hmin*0.1_wp
              realk=.true.
           endif
           call daprop_lr(nch,ninit,nopen,yk,x)
           nldr=nldr+1
           ii=0
           do j=1,nch
              do i=1,j
                 ii=ii+1
                 y(ii)=yk(i,j)
              enddo
           enddo
           safety=.true.   
        endif
        nend=0
        h=hmin
     else
        safety=.true.
        h=hnext
     endif

     ! 007 - handling of possible local resonances 
  enddo cycle_over_r
  return
end subroutine varp_rk
      

subroutine rkqs(y,dydx,n,x,htry,yscal,hdid,hnext,nch,nopen,ninit)
  use precision
  use var_ph_module
  use shared
  implicit none
  integer, intent(in) ::  n, ninit, nch, nopen
  real(kind=wp), dimension(n), intent(inout) :: y, dydx, yscal
  real(kind=wp), intent(in)  :: htry
  real(kind=wp), intent(out) :: hnext, hdid, x

  real(kind=wp) :: h, htemp, xnew, errmax  
  real(kind=wp), dimension(n) :: yerr, ytemp 
  real(kind=wp), parameter :: safety=0.9, pgrow=-.2 
  real(kind=wp), parameter :: pshrnk=-.25
  real(kind=wp), parameter :: errcon=1.89e-4 
  real(kind=wp), parameter :: hcut=2.0
  integer i
  
  !     uses derivs,rkck
  !     fifth-order runge-kutta step with monitoring of local truncation 
  !     error to ensure accuracy and adjust stepsize. input are the 
  !     dependent variable vector y(1:n) and its derivative dydx(1:n) 
  !     at the starting value of the independent variable x. 
  !     also input are the stepsize to be attempted htry, 
  !     the required accuracy eps, and the vector yscal(1:n) against 
  !     which the error is scaled. on output, y and x are replaced 
  !     by their new values, hdid is the stepsize that was actually 
  !     accomplished, and hnext is the estimated next stepsize. 
  !     derivs is the user-supplied subroutine that computes the 
  !     right-hand side derivatives. 
  
  
  !     the value errcon equals (5/safety)**(1/pgrow), 
  !     see use below. 

  h=htry 
1 call rkck(y,dydx,n,x,h,ytemp,yerr,nch,nopen,ninit)
  
  errmax=0. 
  do i=1,n
     errmax=max(errmax,abs(yerr(i)/yscal(i))) 
  enddo
  errmax=errmax/eps 
  if(errmax.gt.1.)then 
     htemp=safety*h*(errmax**pshrnk)
     h=sign(max(abs(htemp),0.1*abs(h)),h) 
     xnew=x+h 
     if(xnew.eq.x) then
        write(*,*) 'stepsize underflow in rkqs'
        stop
     endif
     goto 1 
  else 
     if(errmax.gt.errcon)then
        hnext=safety*h*(errmax**pgrow) 
     else 
        hnext=5.*h 
     endif
      if(care.and.hnext>=hcut)hnext=hcut
     hdid=h 
     x=x+h 
     do  i=1,n
        y(i)=ytemp(i) 
     enddo
     return 
  endif
end subroutine rkqs



subroutine rkck(y,dydx,n,x,h,yout,yerr,nch,nopen,ninit)
  use precision
  use var_ph_module
  implicit none
  integer, intent(in):: n, nch, nopen, ninit
  real(kind=wp), intent(inout), dimension(n) :: dydx, y, yerr, yout 
  real(kind=wp), intent(inout) :: h, x

  
  !     given values for n variables y and their derivatives 
  !     dydx known at x, use the fifth-order cash-karp runge-kutta 
  !     method to advance the solution over an interval h and return 
  !     the incremented variables as yout. also return an estimate 
  !     of the local truncation error in yout using the embedded 
  !     fourth-order method. the user supplies the subroutine 
  !     derivs(x,y,dydx), which returns derivatives dydx at x. 
  integer i   
  real(kind=wp), dimension(n) :: ak2(n),ak3(n),ak4(n),ak5(n),ak6(n),ytemp(n)
  real(kind=wp), parameter :: a2=0.2_wp, a3=0.3_wp, &
       a4=0.6_wp, a5=1.0_wp, a6=0.875_wp, &
       b21=0.2_wp, b31=3.0_wp/40.0_wp, b32=9.0_wp/40.0_wp, b41=0.3_wp, &
       b42=-0.9_wp, b43=1.2_wp, b51=-11.0_wp/54.0_wp, b52=2.50_wp, & 
       b53=-70.0_wp/27.0_wp, b54=35.0_wp/27.0_wp, & 
       b61=1631.0_wp/55296.0_wp, &
       b62=175.0_wp/512.0_wp, b63=575.0_wp/13824.0_wp, &
       b64=44275.0_wp/110592.0_wp, &
       b65=253.0_wp/4096.0_wp, c1=37.0_wp/378.0_wp, c3=250.0_wp/621.0_wp, & 
       c4=125.0_wp/594.0_wp, c6=512.0_wp/1771.0_wp
  real(kind=wp), parameter :: dc1=c1-2825.0_wp/27648.0_wp, &
       dc3=c3-18575.0_wp/48384.0_wp, dc4=c4-13525.0_wp/55296.0_wp, &
       dc5=-277.0_wp/14336.0_wp, dc6=c6-0.25_wp


  do  i=1,n 
     ytemp(i)=y(i)+b21*h*dydx(i) 
  enddo
  call kder(x+a2*h,ytemp,ak2,n,nch,nopen,ninit)
  do  i=1,n
     ytemp(i)=y(i)+h*(b31*dydx(i)+b32*ak2(i)) 
  enddo
  call kder(x+a3*h,ytemp,ak3,n,nch,nopen,ninit)
  do  i=1,n
     ytemp(i)=y(i)+h*(b41*dydx(i)+b42*ak2(i)+b43*ak3(i))
  enddo
  call kder(x+a4*h,ytemp,ak4,n,nch,nopen,ninit)
  do  i=1,n
     ytemp(i)=y(i)+h*(b51*dydx(i)+b52*ak2(i)+b53*ak3(i)+ &
          b54*ak4(i))
  enddo
  call  kder(x+a5*h,ytemp,ak5,n,nch,nopen,ninit)
  do  i=1,n
     ytemp(i)=y(i)+h*(b61*dydx(i)+b62*ak2(i)+b63*ak3(i)+ &
          b64*ak4(i)+b65*ak5(i))
  enddo
  call kder(x+a6*h,ytemp,ak6,n,nch,nopen,ninit)
  do  i=1,n 
     yout(i)=y(i)+h*(c1*dydx(i)+c3*ak3(i)+c4*ak4(i)+ &
          c6*ak6(i))
  enddo
  do  i=1,n
     yerr(i)=h*(dc1*dydx(i)+dc3*ak3(i)+dc4*ak4(i)+dc5*ak5(i) &
          +dc6*ak6(i))
  enddo
  return
end subroutine rkck















