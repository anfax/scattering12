!!!*************************************************************
! 文件/File: close_coupling.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: close_coupling.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

subroutine close_coupling(jtot, step_par, asympt_channels, open_channels, &
     local_channels, symmetry, xs2, myrank, nproc, first_time, &
     multiplicity)
     
  use precision
  ! this subroutine handles directly the scattering event
  ! "local channels" refers to the v,j,l,r labelling   
  ! "channels" refer to the asymptotic channels 
  !
  use shared
  use molecular_parameters
  use potential, only : num_lambda_terms
  use long_range_mod
  use percival_seaton
  use channels

  implicit none
  ! +-------+   
  ! | input |
  ! +-------+
  integer, intent(IN) :: asympt_channels   ! number of asymtotic channels
  integer, intent(IN) :: jtot              ! partial wave index
  integer, intent(IN) :: step_par          ! step parameter
  integer, intent(IN) :: symmetry          ! symmetry identifier (omega)
  integer, intent(IN) :: myrank, nproc     ! mpi variables
  integer, intent(INOUT) :: local_channels ! dimension of the            
                                           ! scattering matrix
                                           !   setted to zero on return 
                                           !   means that
                                           !   we found no open 
                                           !   local channels in
                                           !   the symmetry considered  
  logical, intent(IN) :: first_time        ! first call flag (not used)
  integer, intent(IN) :: multiplicity      ! diatom multiplicity 
                                           ! (s->1, d->2, t->3)

  ! +--------+   
  ! | output |
  ! +--------+
  ! symmetry dependent cross-sections
  ! assumed shape array: see interface body in main()
  real(kind=wp), dimension(:,:,:,:,:,:), intent(OUT) :: xs2
  ! number of asymptotic open channels |n,j,r> 
  integer ::  open_channels 
  
  ! declaration of local variables follows

  ! number of open local channels (n,j,l,r) 
  integer       :: open_local_channels, ich
  ! functions for 3-j and 6-j coefficients
  real(kind=wp) :: xf3j, xf6j
  ! Log-Derivative Matrix
  real(kind=wp), dimension(local_channels,local_channels) :: y
  ! T-matrix
  real(kind=wp), allocatable, dimension(:,:) :: ti, tr, t_square 

  integer       :: counter, i, ii, l, j1, j2, l1, l2, n1, n2 
  integer       :: j, r1, r2, unit 
  integer       :: s1, s2
  integer       :: lmin, lmax, lambda_step, lambda_term, lambda 
  real(kind=wp) :: zero=0.0_wp, one=1.0_wp, two=2.0_wp, sign, coeff, over
  ! the following are all variables used in the eventual change of rend
  logical       :: rend_changed  
  real(kind=wp) :: new_rend, desired_centrifugal, desired_rend, vlast
  real(kind=wp) :: max_wavevector

  !unit on which we write debug information
  unit=myrank+20   
  if (nproc==1) unit=6  ! std output

  ! we do not need to count channels, but we need to identify 
  ! the two parity cases and assign n_qn j_qn and lorb vectors.
  if (debug>=1) then
     write(unit,1041)'** Close_Coupling: called with jtot:',jtot,'n_ch:', &
          local_channels,'symm:',symmetry
  endif
1041 format(t20,a,2x,i4,5x,a,2x,i5,5x,a,2x,i3)

  ! allocate percival and seaton coeff.
  allocate(pot_coeff(num_lambda_terms,local_channels,local_channels))

  ! allocate quantum numbers and energies
  allocate(     j_local(local_channels))
  allocate(     n_local(local_channels))
  allocate(     l_local(local_channels))
  allocate(     r_local(local_channels))
  allocate(     s_local(local_channels))
  allocate(local_energy(local_channels))
  allocate(wave_vector (local_channels))
    
  local_energy=0.0_wp
  j_local=0.0_wp
  n_local=0.0_wp
  l_local=0.0_wp
  r_local=0.0_wp
  s_local=0.0_wp
  pot_coeff=0.0_wp

  ! start labelling the quantum states
  ! CC equations  
  counter=0
  do i=1,asympt_channels
        lmin=abs(jtot-j_qn(i))
        lmax=jtot+j_qn(i)
        l_cycle: do l=lmin,lmax
           if ((symmetry==1).and.(.not.parity( l+r_qn(i) ))) then 
              counter=counter+1
              n_local(counter)= n_qn(i)
              j_local(counter)= j_qn(i)
              r_local(counter)= r_qn(i) 
              s_local(counter)= s_qn(i)
              l_local(counter)= l
              local_energy(counter)= e_level(i)
           endif
           if ((symmetry==2).and.(parity( l+r_qn(i) ))) then 
              counter=counter+1
              n_local(counter)= n_qn(i)
              j_local(counter)= j_qn(i)
              r_local(counter)= r_qn(i) 
              s_local(counter)= s_qn(i)
              l_local(counter)= l
              local_energy(counter)= e_level(i)
           endif
        enddo l_cycle
  enddo

  ! code below is a check: counter *must* be equal to n_local_channels
  if (counter/=local_channels) then
     write(*,*) 'error in close_coupling: severe inconsistency:'
     write(*,*) 'number of local channels(v,j,l,r) in close_coupling()', counter 
     write(*,*) 'is different from that calculated in main()',local_channels
     stop
  endif

  ! channels output
  if (debug>=2)then
     write(unit,1042)'    writing out channels for:','jtot','msym'
     write(unit,fmt='(t50,i4,t60,i3)')jtot,symmetry
     write(unit,1043)'index','n-vib','j-rot','l-orb','N-rot','energy'  
  endif
  if (debug>=2)then
     do i=1,local_channels
        write(unit,1044)i,n_local(i),j_local(i),l_local(i),r_local(i),local_energy(i)/cm_to_au
     enddo
  endif

  ! now calculate the potential matrix coefficient
  ! note: lambda_term goes on the symmetries of the v_lambda.
  ! vmatrix=sum_over_lambda (v_lambda(r)*c(j,l,r,j',l',r',lambda))
  ! we want the c coefficient
  if (debug>=1)then
     write(unit,999)'** Close_Coupling: calculating v-coeff'
  endif
  if (debug>=3)then
     write(unit,999)'close_coupling: writing out coeff'
     write(unit,1043)'i','ii','lambda','vcoef'
  endif
  
  lambda_step=1
  if (homonuclear) lambda_step=2
  outer_local_channels: do i=1,local_channels
     inner_local_channels: do ii=i,local_channels
        loop_lambda_terms: do lambda_term=1,num_lambda_terms
           lambda=lambda_step*(lambda_term-1)

            ! CC equations
            three_cases: if(multiplicity==1)then
            ! singlet
            ! taken from A. E. DePristo and M. H. Alexander,
            ! J. Phys. B: Atom. Molec. Phys. 9, L39-L41 (1976);
            ! page L40, formula (3).

              j1=j_local(i)
              j2=j_local(ii)
              l1=l_local(i)
              l2=l_local(ii)
              sign=(-1)**(j1+j2-jtot)
              coeff=real((2*j1+1)*(2*j2+1)*(2*l1+1)*(2*l2+1),kind=wp)
              coeff=sign*dsqrt( coeff ) * &
              xf3j(real(j1,kind=wp),real(lambda,kind=wp),real(j2,kind=wp), &
                   zero,            zero,                zero           )* &
              xf3j(real(l1,kind=wp),real(lambda,kind=wp),real(l2,kind=wp), &
                   zero,            zero,                zero           )* &
              xf6j(real(j1,kind=wp),real(j2,kind=wp),real(lambda,kind=wp), &
                   real(l2,kind=wp),real(l1,kind=wp),real(jtot,kind=wp) )
              !--------------------------------
              pot_coeff(lambda_term,i,ii)=coeff
              pot_coeff(lambda_term,ii,i)=coeff
              if ((debug>=3).and.(i<=ii)) then
                 write(unit,1045)i,ii,lambda,coeff
              endif
  
            elseif(multiplicity==2)then
            ! doublet
            ! taken from G. C. Corey and F. C. McCourt,
            ! J. Phys. Chem. 87, 2723-2730 (1983);
            ! page 2725, formula (2.22).           
  
              j1=j_local(i)
              j2=j_local(ii)
              l1=l_local(i)
              l2=l_local(ii)
              r1=r_local(i)
              r2=r_local(ii)
              ! sign=(-1)**(S-LAMBDA-JTOT) 
              sign=(-1)**(1-lambda-jtot) ! 1=S+1/2; JTOT=J(VERO)+1/2
              coeff=(two*j1)*(two*j2)*(two*l1+one)*(two*l2+one)* &
                    (two*r1+one)*(two*r2+one) ! J1=J1(VERO)+1/2, SAME FOR J2
              coeff=sign*dsqrt( coeff ) * &
                   xf3j(real(r2,kind=wp),real(lambda,kind=wp),real(r1,kind=wp), &
                        zero,            zero,                zero           )* &
                   xf3j(real(l2,kind=wp),real(lambda,kind=wp),real(l1,kind=wp), &
                        zero,            zero,                zero           )* & 
                   xf6j(real(j1,kind=wp)-.5d0,real(j2,kind=wp)-.5d0,real(lambda,kind=wp), & 
                        real(l2,kind=wp),real(l1,kind=wp),real(jtot,kind=wp)-.5d0 )* & 
                   xf6j(real(j1,kind=wp)-.5d0,real(j2,kind=wp)-.5d0,real(lambda,kind=wp), &  
                        real(r2,kind=wp),real(r1,kind=wp),        .5d0        )
              !--------------------------------
              pot_coeff(lambda_term,i,ii)=coeff
              pot_coeff(lambda_term,ii,i)=coeff
              if ((debug>=3).and.(i<=ii)) then
                 write(unit,1045)i,ii,lambda,coeff
              endif

            elseif(multiplicity==3) then
            ! triplet
            ! taken from G. C. Corey and F. C. McCourt,
            ! J. Phys. Chem. 87, 2723-2730 (1983);
            ! page 2725, formula (2.22).

              j1=j_local(i)
              j2=j_local(ii)
              l1=l_local(i)
              l2=l_local(ii)
              r1=r_local(i)
              r2=r_local(ii)
              sign=(-1)**(1-lambda-jtot)
              coeff=(two*j1+one)*(two*j2+one)*(two*l1+one)*(two*l2+one)* &
                    (two*r1+one)*(two*r2+one)
              coeff=sign*dsqrt( coeff ) * &
                   xf3j(real(r2,kind=wp),real(lambda,kind=wp),real(r1,kind=wp), &
                        zero,            zero,                zero           )* &
                   xf3j(real(l2,kind=wp),real(lambda,kind=wp),real(l1,kind=wp), &
                        zero,            zero,                zero           )* &
                   xf6j(real(j1,kind=wp),real(j2,kind=wp),real(lambda,kind=wp), &
                        real(l2,kind=wp),real(l1,kind=wp),real(jtot,kind=wp) )* & 
                   xf6j(real(j1,kind=wp),real(j2,kind=wp),real(lambda,kind=wp), &
                        real(r2,kind=wp),real(r1,kind=wp),        one        )
              !------------------------------------------------------------------
              pot_coeff(lambda_term,i,ii)=coeff
              pot_coeff(lambda_term,ii,i)=coeff
              if ((debug>=3).and.(i<=ii)) then
                 write(unit,1045)i,ii,lambda,coeff
              endif

            endif three_cases

              if ((debug>=3).and.(i<=ii)) then
                 write(unit,1045)i,ii,lambda,coeff
              endif 

        enddo loop_lambda_terms
     enddo inner_local_channels
  enddo outer_local_channels
  if (debug>=1)then
     write(unit,999)'** Close_Coupling: ok, we have the v-coeff matrix'
  endif

  !       if we see that some channels are under the centrifugal
  !       barrier at large distances for some jtot, if the last 
  !       point of the potential doesn't satisfy this condition
  !       program stops.
  open_channels=0
  open_local_channels=0
  ! this is for the open thriples (j,l,n,r)
  do i=1,local_channels
     if (local_energy(i)<=total_energy) then
        open_local_channels=open_local_channels+1
     endif
  enddo

  if (debug>=1) then
     write(unit,fmt='(t20,a,2xi4)')'** Open (n,j,l,r) channels are:',open_local_channels
  endif

  ! this is for the open couples (j,n,r)
  do i=1,asympt_channels
     if (e_level(i)<=total_energy) then
        open_channels=open_channels+1
     endif
  enddo
  if (debug>=1) then
     write(unit,fmt='(t20,a,2xi4)')'** Open (n,j,r)  channels are:',open_channels
  endif

  ! if there are no open_local_channels return with no cross sections and
  ! set the local_channel flag for correct handling of the situation in main()
  if (open_local_channels==0) then
     local_channels=0
     if (debug>=1) then
        write(*,fmt='(t20,a,i4)')'No open (local) channels for symmetry ', symmetry
        write(*,fmt='(t20,a)')'passing control to main() with no propagation'
     endif
     deallocate(pot_coeff)
     deallocate(j_local,n_local,l_local,r_local,local_energy, wave_vector)
     return
  endif

  ! check if the centrifugal barrier masquerade some of the open channel
  if (.not.use_variable_phase) rswitch=rend 

  over=1.0/((rend*rend)*two_mu)
  new_rend=rend
  rend_changed=.false.
  do i=1,open_local_channels
     vlast=l_local(i)*(l_local(i)+1)*over

     if ((total_energy-local_energy(i))<vlast) then
        if (irxset<=0.0_wp) then
           write(*,*)'Close_coupling: warning!'
           write(*,*)'    channel energy(cm-1):',local_energy(i)/cm_to_au
           write(*,*)'                         ',(total_energy-local_energy(i))/cm_to_au
           write(*,*)'    cent. potential(cm-1):',vlast/cm_to_au
           write(*,*)'  there are some channels under the centrifugal'
           write(*,*)'  barrier at long range: increase rend'
           write(*,*)'  program stopped'
           stop
        else
           rend_changed=.true.
           desired_centrifugal=real(l_local(i)*(l_local(i)+1),kind=wp) &
                /(two_mu*(total_energy-local_energy(i)))
           desired_rend=irxset*sqrt(desired_centrifugal)
           if (desired_rend>new_rend) then 
              new_rend=desired_rend
              ich=i
           endif
        endif
     endif
  enddo
  if(rend_changed) then
     write(unit,999) 'Close_coupling: warning!'
     write(unit,*)'                   ','channel: ',ich,n_local(ich),j_local(ich),l_local(ich),r_local(ich)
     write(unit,1046)'    channel energy(cm-1) :',local_energy(ich)/cm_to_au
     write(unit,1046)'    etot-ch_energy(cm-1) :',(total_energy-local_energy(ich))/cm_to_au
     write(unit,1046)'    cent. potential(cm-1):',vlast/cm_to_au
     write(unit,999) 'some of the open channels are below the centrifugal barrier'
     rend=new_rend
     write(unit,1046)'at rend, REND INCREASED TO: ',rend*au_to_angst
     if (.not.use_variable_phase) then
        rswitch=new_rend
     endif
  endif

  ! calculating number of steps for log-derivative part of propagation
  ! nsteps=(rswitch-rbegin) / [pi/k_max*steps)]
  ! we do this only for the first energy 
  !if(first_time) then
  !      max_wavevector= sqrt( two_mu * (total_energy-e_level(1)) )
  !      step_size=pi/(max_wavevector*step_par)
  !      nsteps=(rswitch-rbegin)/step_size
        nsteps=step_par
        step_size=(rswitch-rbegin)/nsteps
        half_step_size=step_size/2.0_wp
  !endif

  if (debug>=1) then
     write(unit,fmt='(t20,a,1x,i6)')'** Number of step for pure Log-derivative' &
          ,nsteps
     write(unit,fmt='(t20,a,f15.8)')'** Stepsize: ', step_size
  endif

1042 format(t20,a,t50,a,t60,a)
1043 format(t20,a,t30,a,t40,a,t50,a,t60,a,t70,a)
1044 format(t20,i3,t30,i3,t40,i3,t50,i3,t60,i3,t70,f8.2)
1045 format(t20,i4,t30,i4,t40,i4,t50,f8.4)
1046 format(t20,a,2x,f16.8,5x,a,i6)

 ! allocate and initialize T-matrix
  allocate(tr      (open_local_channels,open_local_channels))
  allocate(ti      (open_local_channels,open_local_channels))
  allocate(t_square(open_local_channels,open_local_channels))
  ti=0.0_wp
  tr=0.0_wp
  t_square=0.0_wp

  !   *************************************
  !   +-----------------------------------+
  !   |       Calling Propagators         |
  !   +-----------------------------------+
  !   *************************************

  write(unit,fmt= '(a,i3,t10,a,a,f10.3,a,f10.3,a,i7,a)') &
       'RK-',myrank, &
       '-> Log-Derivative propagation: ', &
       'from ',rbegin*au_to_angst,'  to ',rswitch*au_to_angst, &
       ' in ',nsteps, ' steps'

  ! multiply energies by 2*mu (log-der wants k-vectors)  
  local_energy=local_energy*two_mu
  wave_vector=sqrt(abs(reduced_e-local_energy))
  if(use_variable_phase)then

   call logdvp(local_channels, open_local_channels, tr, ti, myrank, nproc)

  else

   call logd(local_channels, open_local_channels, tr, ti, myrank, nproc)

  endif

  ! in y there is the log-derivative Y matrix at rswitch
  
  !restoring energies for next cycle
  local_energy=local_energy/two_mu

  t_square(1:open_local_channels,1:open_local_channels)=     &
       +ti(1:open_local_channels,1:open_local_channels)**2   &
       +tr(1:open_local_channels,1:open_local_channels)**2 
  ! sum over all the thriple j,n,r =j',n',r' the different l values

  do i=1,open_local_channels
     j1=j_local(i)+1
     n1=n_local(i)+1
     r1=r_local(i)+2  
     s1=s_local(i)
     do ii=1,open_local_channels
        j2=j_local(ii)+1
        n2=n_local(ii)+1
        r2=r_local(ii)+2
        s2=s_local(ii)

        xs2(n1,n2,j1,j2,s1,s2)=xs2(n1,n2,j1,j2,s1,s2)+ &
             ( t_square(i,ii)/real((2*j_local(i) + 1),kind=wp) )

     enddo
  enddo
999 format(t20,a)
  deallocate(tr,ti,t_square)
  deallocate(pot_coeff)
   deallocate(j_local,n_local,l_local,r_local,s_local,local_energy, &
              wave_vector)
  return
end subroutine close_coupling

