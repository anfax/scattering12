!!!*************************************************************
! 文件/File: kmatrix.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: kmatrix.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


subroutine kmatrix(J,Jp,N,Np,charge,energy)
  use Solve,only : Z=>eigen_vector,b=>eigen_value,invert_real,invert_complex,diag_complex,diag_symm
  use Open_information,only : max_open,kappa,E,max_index,lmax,wf_ind,max_open_total
  use Brep,only : R0,p_Patch,x_coord=>coord_x,y_coord=>coord_y,z_coord=>coord_z,boundary,plotpoints
  use nrtype, only :Pi,dbl,i4b,dpc
  use bsp_info
  use control , only : output_file,molecule,option_wf
  use Matrices, only : RHS_co => d
  use V_setup, only : n1,n2,n3
  implicit none
  interface perturbation_87
     subroutine perturbation_87(max_open,R0,charge,energy,lmax,K_perturbed,Z_mat,b)
       use nrtype,only : i4b,dbl,dpc,Pi
       integer(kind=i4b),intent(in) :: max_open,lmax
       real(kind=dbl),intent(in) :: charge,energy,R0,Z_mat(max_open,max_open),b(max_open)
       real(kind=dbl),intent(out) :: K_perturbed(max_open,max_open)
     end subroutine perturbation_87
  end interface
  integer(kind=i4b) :: p,q,f,indice,alpha,enne,nClosed,max_9000
  integer(kind=i4b),dimension(:),allocatable :: indexing
  integer(kind=i4b),save :: counter=0
  integer(kind=i4b),parameter :: nx=13,ny=14,nz=15
  real(kind=dbl),dimension(max_open) :: N4,J,Jp,N,Np,tan_qd
  real(kind=dbl),dimension(max_open) ::sigma_partial
  real(kind=dbl)::sigma_total,sigma_elastic,sigma_inelastic,sigma_elastic2,sigma_elastic3,sigma1,sigma2,Ene,eigenphase_sum,delta,eta,eigenphase_sum_2,charge,energy
  real(kind=dbl),dimension(:,:),allocatable :: R,K,M,T,mostro,astro,identity,I_mat,J_mat,F_mat,Fp_mat,Z_mat,cc,U_vec,K_perturbed
  real(kind=dbl),dimension(:),allocatable :: W
  real(kind=dbl),parameter :: angstrom_atomicu_conversion=0.5291772083d0
  real(kind=dbl) :: x,y,zzz
  complex(kind=dpc),dimension(:,:),allocatable :: L,S,ident,probe,cc_complex,cc_complex2,grid_func
  complex(kind=dpc),dimension(:,:),allocatable,save :: coeff_projection
  complex(kind=dpc),dimension(:),allocatable :: diag_S,grid
  complex(kind=dpc),parameter :: I=(0.0d0,1.0d0) 
  complex(kind=dpc) :: det,det2,argom,gamma_s,gammaln,cgamma
  character*10 :: method='No perturbation'
  external :: gammaln,cgamma

  !---------------------------
  ! Subroutine that calculates the scattering K,S matrices, cross sections and
  ! quantum defects
  !---------------------------

  allocate(R(max_open,max_open))
  allocate(K(max_open,max_open))
  allocate(M(max_open,max_open))
  allocate(L(max_open,max_open))
  allocate(S(max_open,max_open))

  allocate(Z_mat(max_open,max_open))

  !**COMM**********
  !Wavefunction plotting



  nClosed=p_Patch-max_index
  write(6,*)'nClosed=',nClosed,'p_Patch=',p_Patch,'max_index=',max_index
  !  write(66,*)nClosed,p_Patch,max_index
  write(66,*)plotpoints*8+max_index,plotpoints*8+2*max_index,max_index
  allocate(cc(nClosed,max_open))
  ! Debug
  !  cc= - matmul(RHS_co,Z)

  ! Call Blas 3 subroutine to perform this multiplication
  cc(:,:)=0.d0
  call dgemm ( 'N', 'N', p_Patch-max_index, max_open, max_open, -1.d0, RHS_co, p_Patch-max_index, Z, max_open, &
       &                   0.d0, cc, p_Patch-max_index )

  !cc=cc*sqrt(2.d0)
  ! Wf plotting
  !allocate(indexing(nClosed))
  !q=0
  !do p=1,nClosed
  !if (mod(wf_ind(p)-1,8)==0) then
  !q=q+1
  !indexing(q)=p
  !end if
  !end do
  !max_9000=q
  !--------------------
if (option_wf=='wf') then
! Write the wavefunctions (R-matrix normalized) to file
  f=1
  q=0
  do f=1,max_open
     q=0
     !     do p=1,nClosed-max_index,8
     do p=1,plotpoints*8,8
        q=q+1
        !do q=1,max_9000
                write(9000+f,902)x_coord(q),y_coord(q),z_coord(q),cc(p,f)
        !       write(9000+f,902)x_coord(q),y_coord(q),z_coord(q),cc(indexing(q),f)
     end do
  end do
end if
902 format(5e13.5)
  deallocate(cc)
  !deallocate(indexing)
  !**COMM**********

  N4=0.0d0
  R=0.0d0
  K=0.0d0
  M=0.0d0
  L=cmplx(0.0d0,0.0d0)
  S=cmplx(0.0d0,0.0d0)
  kappa=sqrt(2*E)

  do p=1,max_open
     do q=1,max_open
        do f=1,max_open
           N4(p)=N4(p)+Z(f,q)*Z(f,p)
           M(p,q)=M(p,q)+Z(f,q)*Z(f,p)
        end do
     end do
     if (N4(p).lt.0) then
        N4(p) = -N4(p)
     end if
  end do


  ! R-matrix calculation
  write(6,*)'Rmatrix'


  R=0.0d0
  allocate(mostro(max_open,max_open))
  allocate(T(max_open,max_open))
  !  write(6,*)'after alloc',max_open
  T=0.0d0
  mostro=0.0d0
  do p=1,max_open
     do q=1,max_open
        do f=1,max_open
           R(p,q)=R(p,q)+Z(p,f)*Z(q,f)/b(f)
        end do
     end do
  end do
  !  write(6,*)'after /'
  do p=1,max_open
     do q=1,max_open
        write(63,*)Z(p,q),p,q
        write(76,*)R(p,q),p,q
     end do
     write(64,*)b(p),p
  end do

  ! Matching with asymptotic functions
  ! Invert (N-N'R)

  T=0.0d0
  mostro=0.0d0

  do p=1,max_open
     do q=1,max_open
        if (p.eq.q) then
           T(p,q) = N(p)
        else
           T(p,q) = 0.d0  
        end if
     end do
  end do

  !write(6,*)'1st inversion'
  Z_mat=Z
  Z=0.0d0

  do p=1,max_open
     do q=1,max_open
        M(p,q) = T(p,q) - Np(p)*(R(p,q))     
     end do
  end do
  call write_simm(M,max_open,51)

  T=0.0d0
  mostro=0.0d0
  T=M
  !write(6,*)'invert_real'
  call invert_real(M,max_open,max_open)
  call write_simm(M,max_open,53)

  mostro=matmul(T,M)

  ! Calculate  J-J'R

  T=0.0d0
  mostro=0.0d0

  do p=1,max_open
     do q=1,max_open
        if (p.eq.q) then
           T(p,q) = J(p)
           mostro(p,q)=J(p)/Jp(p)
        else
           T(p,q) = 0
        end if
     end do
  end do

  Z=0.0d0

  do p=1,max_open
     do q=1,max_open
        Z(p,q) = T(p,q) - Jp(p)*R(p,q)
     end do
  end do
  call write_simm(Z,max_open,52)

  ! K-matrix 

  write(6,*)'Kmatrix'

  K=matmul(Z,M)
  do p=1,max_open
     do q=1,max_open
        write(599,*)K(p,q),p,q
     end do
  end do
  !*COMM*****
  !*Diagonalize K-matrix to obtain the quantum defects
  !write(6,*)'diag K'
  T=K


  call diag_symm(T,max_open,max_open,tan_qd)

  do f=1,max_open
     write(103,*)(kappa**2)/2,f,tan_qd(f)
  end do
  !Calculate the quantum defects

  write(6,*)'QD/phase shift modulo 1 at E=',(kappa**2)/2

  if (counter.eq.0)  then
    open(unit=nx,file='phase_shifts.out',status='unknown',action='write')
  else
    open(unit=nx,file='phase_shifts.out',status='old',action='write',access='append')
  end if
  do f=1,max_open
     tan_qd(f)=atan(tan_qd(f))/Pi
     write(474,*)(kappa**2)/2,tan_qd(f),f
     write(nx,*)(kappa**2)/2,tan_qd(f),f
  end do
  close(nx)
  write(6,*)tan_qd(:)

  ! Calculate the eigenphase sums
  do f=1,max_open
    if (tan_qd(f).lt.0.d0) tan_qd(f)=tan_qd(f)+1.d0
    eigenphase_sum=eigenphase_sum+tan_qd(f)
  end do
  write(6,*)'eigenphase sum=',(kappa**2)/2.d0,eigenphase_sum
  write(475,*)(kappa**2)/2.d0,eigenphase_sum



  !-------------------------------------------------------------------------------------
  ! Perturbative (multipole) part (does not work in this version of the code)
  if (method=='Sea99') then
     allocate(K_perturbed(max_open,max_open))
     call perturbation_99(max_open,K,K_perturbed,charge,energy,R0)
     deallocate(K_perturbed)
  else if (method=='Sea87') then
     allocate(K_perturbed(max_open,max_open))
     call perturbation_87(max_open,R0,charge,energy,lmax,K_perturbed,Z_mat,b)
     call get_qd(K_perturbed,max_open,tan_qd,kappa)
     deallocate(K_perturbed)
  else
     write(6,*)'No perturbative treatment'
  end if
  !-------------------------------------------------------------------------------------


  ! Write eigenvectors of K-matrix to file
  write(472,*)max_open
  do p=1,max_open
     do q=1,max_open
        write(472,*)T(p,q)
     end do
  end do

  ! Calculate S-matrix
  L=cmplx(0.0d0,0.0d0)
  allocate(probe(max_open,max_open))
  allocate(identity(max_open,max_open))
  identity=0.0d0
  do p=1,max_open
     do q=1,max_open
        if (p.eq.q) then
           identity(p,q)=1
        else
           identity(p,q)=0
        end if
     end do
  end do

  do p=1,max_open
     do q=1,max_open
        L(p,q) = cmplx(identity(p,q),-K(p,q))
     end do
  end do

  allocate(ident(max_open,max_open))
  probe=cmplx(0.0d0,0.0d0)
  probe=L
  ident=cmplx(0.0d0,0.0d0)

  call invert_complex(L,max_open,max_open)
  !write(6,*)'after invert complex'
  ident=matmul(probe,L)

  S=cmplx(0.0d0,0.0d0)
  S=matmul(cmplx(identity,K),L)
  do p=1,max_open
     do q=1,max_open 
        write(537,*)S(p,q),p,q
     end do
  end do


  !*******************************************************
  ! Calculate I,J matrices
  Z=Z_mat
  allocate(I_mat(max_open,max_open))
  allocate(J_mat(max_open,max_open))
  allocate(W(max_open))
  do p=1,max_open
     W(p)=1/(N(p)*Jp(p)-Np(p)*J(p))
     write(601,*)W(p)
  end do
  do p=1,max_open
     do q=1,max_open
        I_mat(p,q) = (N(p)*Z(p,q)*b(q)-Np(p)*Z(p,q)) * (W(p))
        J_mat(p,q) = (J(p)*Z(p,q)*b(q)-Jp(p)*Z(p,q)) * (W(p))
     end do
  end do
  M=I_mat

!------------------------------------------------------------------------------------------------------------------
! Photoionization Code
if (molecule.eq.'ion') then
  allocate(cc(nClosed,max_open))
!!$ Transform to F solution (eq. 2.30 orange review)
  !  cc= - matmul(RHS_co,Z)

  ! Call Blas 3 subroutine to perform this multiplication
  cc(:,:)=0.d0
  call dgemm ( 'N', 'N', p_Patch-max_index, max_open, max_open, -1.d0, RHS_co, p_Patch-max_index, Z, max_open, &
       &                   0.d0, cc, p_Patch-max_index )
  allocate(cc_complex2(nClosed,max_open))

!!$ Convert to S-matrix eigenstates (eq. 2.41 orange review), with outgoing waves  boundary conditions
  ! cc_complex2=matmul(cc,inverse(I_mat-I*J_mat))
  cc_complex2(:,:)=cmplx(cc(:,:),0.d0) 
  deallocate(cc)
  allocate(cc_complex(nClosed,max_open))
!--------------------------------------------------------------------
!!$ Outgoing BC -----  photoionization
!!$  probe=I_mat+I*J_mat
!!$ Incoming BC -----  HHG experiment
   probe=I_mat-I*J_mat
!--------------------------------------------------------------------

  call invert_complex(probe,max_open,max_open)

  call zgemm ( 'N', 'N', p_Patch-max_index, max_open, max_open, (1.d0,0.d0), cc_complex2, p_Patch-max_index, probe, max_open, &
       &                   (0.d0,0.d0), cc_complex, p_Patch-max_index )
deallocate(cc_complex2)
  f=1
  q=0
  do f=1,max_open
     q=0
     do p=1,plotpoints*8,8
        q=q+1
        !write(9000+f,902)x_coord(q),y_coord(q),z_coord(q),real(cc_complex(p,f)),aimag(cc_complex(p,f))
     end do
  end do

  ! Call Psi^- interpolation
  allocate(grid_func(max_open,n1*n2*n3))
  allocate(grid(n1*n2*n3))
  !write(6,*)'before plot'
  call plot_grid1(cc_complex,grid_func,nClosed,p,J,N,S,max_open)
  !write(6,*)'after plot'
end if
  ! Call projection of numerical harmonics
  if (counter.eq.0) then
  write(6,*)'after plot',counter
  if (molecule.eq.'ion')deallocate(cc_complex) ! test a bit
  allocate(coeff_projection(max_open_total,max_open_total))
  call harmonic_projection1(cc_complex(:,p),grid_func,nClosed,p,J,N,S,max_open_total,coeff_projection)
  end if
  counter =counter+1
if (molecule.eq.'ion') then
  ! Call integration and photoionization cross section

  call assemble_routine1(cc_complex(:,p),grid_func,nClosed,p,J,N,S,max_open,coeff_projection(1:max_open,1:max_open),energy)

 ! deallocate(coeff_projection)
  deallocate(grid_func,grid)
end if ! End photoionization code
!--------------------------------------------------------------------------------------------------------------------------
  !T=K
  !call diag_symm(T,max_open,max_open,tan_qd)



  !sigma_elastic=0.d0
  !do f=1,max_open
  !   tan_qd(f)=atan(tan_qd(f))/Pi
  !   write(490,*)(kappa**2)/2,f,tan_qd(f)
  !   sigma_elastic=sigma_elastic+Pi/(kappa*kappa)*4*sin(Pi*tan_qd(f))**2
  !end do
  !write(701,*)sigma_elastic,kappa**2/2


  !*******************************************************
  ! ** Transformations to plot the wavefunction -- According to Chris'
  ! ** suggestions for ions

  !do p=1,max_open
  !   do q=1,max_open
  !      T(p,q)=T(p,q)*cos(Pi*tan_qd(q))
  !   end do
  !end do

  !cc=matmul(cc,I_mat)
  !cc=matmul(cc,T)

  !do f=1,max_open
  !   q=0
  !   do p=1,nClosed-max_index,8
  !      q=q+1
  !      write(6,*)q
  !      write(10000+f,902)x_coord(q),y_coord(q),z_coord(q),cc(p,f)
  !   end do
  !end do
  !*******************************************************
  !Calculate F matrix

  !allocate(F_mat(max_open,max_open))
  !allocate(Fp_mat(max_open,max_open))

  !I_mat=M

  !do p=1,max_open
  !   do q=1,max_open
  !      F_mat(p,q)=J(p)*I_mat(p,q)-N(p)*J_mat(p,q)
  !      Fp_mat(p,q)=Jp(p)*I_mat(p,q)-Np(p)*J_mat(p,q)
  !   end do
  !end do


  !call invert_real(Fp_mat,max_open,max_open)

  !R=matmul(F_mat,Fp_mat)

  !do p=1,max_open
  !   do q=1,max_open
  !      write(64,*)R(p,q),p,q
  !   end do
  !end do


  !T=0.0d0
  !mostro=0.0d0

  !do p=1,max_open
  !   do q=1,max_open
  !      if (p.eq.q) then
  !         T(p,q) = N(p)
  !      else
  !         T(p,q) = 0
  !      end if
  !   end do
  !end do

  !write(6,*)'1st inversion'

  !Z=0.0d0

  !do p=1,max_open
  !   do q=1,max_open
  !      M(p,q) = T(p,q) - Np(p)*(R(p,q))
  !   end do
  !end do
  !T=0.0d0
  !mostro=0.0d0
  !T=M
  !write(6,*)'before invert_real'
  !call invert_real(M,max_open,max_open)
  !write(6,*)'before mostro'
  !mostro=matmul(T,M)
  !T=0.0d0
  !mostro=0.0d0

  !do p=1,max_open
  !   do q=1,max_open
  !      if (p.eq.q) then
  !         T(p,q) = J(p)
  !         mostro(p,q)=J(p)/Jp(p)
  !      else
  !         T(p,q) = 0
  !      end if
  !   end do
  !end do

  !Z=0.0d0

  !do p=1,max_open
  !   do q=1,max_open
  !      Z(p,q) = T(p,q) - Jp(p)*R(p,q)
  !   end do
  !end do

  !write(6,*)'K'



  !*******************************************************

  !Calculate Eigenphase sums

  !eigenphase_sum=0.d0
  !eigenphase_sum_2=0.d0
  !do f=1,max_open
  !   eigenphase_sum=eigenphase_sum+tan_qd(f)
  !   eigenphase_sum_2=eigenphase_sum_2+abs(tan_qd(f))
  !end do
  !write(475,*)eigenphase_sum,(kappa**2)/2
  !write(469,*)eigenphase_sum_2*Pi,(kappa**2)/2

  !*COMM*******
  !*Calculation of the energies

  !do alpha=1,max_open
  !   do enne=1,max_open
  !      Ene=-1.d0/(2.d0*(enne-tan_qd(alpha))**2)
  !      write(473,*)Ene,enne,alpha
  !   end do
  !end do

  !*COMM*****


  ! Diagonalization of the S matrix
  allocate(diag_S(max_open))
!-----------------------------------
! S-matrix projection on complex spherical harmonics
ident=S
ident=matmul(ident,transpose(conjg(coeff_projection/R0)))
ident=matmul(coeff_projection/R0,ident)
  do p=1,max_open
     do q=1,max_open 
        write(536,*)ident(p,q),p,q
     end do
  end do

diag_S=cmplx(0.d0,0.d0)
  call diag_complex(ident,L,max_open,max_open,diag_S)
  !write(6,*)'after diag_complex'
  sigma_elastic=0.0d0
  write(6,*)'cross sections'
  do p=1,max_open
     sigma_elastic=sigma_elastic+Pi/(kappa*kappa)*(diag_S(p)-1)*conjg(diag_S(p)-1) 
  end do
  write(6,*)'cross section xform=',sigma_elastic, kappa**2 / 2
! Verify xform is unitary
!ident=S
!ident=matmul(ident,transpose(conjg(coeff_projection/R0)))
!ident=matmul(coeff_projection/R0,ident)
!write(36,*)'Smataa',matmul(ident,transpose(conjg(ident)))
!-------------------------------------------------------
!-----------------------------------

  !***COMM***
  !*Diagonalize S-matrix

  call diag_complex(S,L,max_open,max_open,diag_S)
!write(6,*)acos(real(diag_S(:)))/Pi !! TEST

  !write(6,*)'after diag_complex'
  !do p=1, max_open
  !   do q=1,max_open
  !      write(399,*)S(p,q),p,q
  !   end do
  !   write(398,*)diag_S,p
  !end do
  !do p=1,max_open
  !   write(39,*)real(1/I*(S(p,p))),p
  !end do

  !***COMM****************************
  !*Calculate the Coulomb phase shift
  !do p=1,max_open
  !   argom=1.d0*p+1.d0+I/kappa
  !   gamma_s=cgamma(argom)
  !   eta=atan(aimag(gamma_s)/ real(gamma_s))
  !   write(104,*)eta,(kappa**2)/2,p
  !end do
  !***COMM****************************
 ! write(6,*)'eigenphase sum'
  eigenphase_sum=0.d0
  do p=1,max_open
     delta=.5d0*acos(1.-real((diag_S(p)-1)*conjg(diag_S(p)-1))/2.d0)
     eigenphase_sum=eigenphase_sum+delta
  end do
  write(478,*)eigenphase_sum,(kappa**2)/2



  sigma_elastic=0.0d0
  !write(6,*)'cross sections'
  ! Elastic
  do p=1,max_open
     sigma_elastic=sigma_elastic+Pi/(kappa*kappa)*(diag_S(p)-1)*conjg(diag_S(p)-1) 
  end do
  write(output_file,*) kappa**2 / 2,sigma_elastic
  write(6,*)'cross section=',sigma_elastic,kappa**2 / 2

  !deallocate(cc_complex)
  deallocate(identity)
  deallocate(probe)
  deallocate(ident)
  deallocate(diag_S)
  deallocate(mostro)
  deallocate(T)
  deallocate(R)
  deallocate(K)
  deallocate(M)
  deallocate(L)
  deallocate(S)
  ! New deallocations
  deallocate(Z_mat)
  !deallocate(I_mat)
  !deallocate(J_mat)
  !deallocate(F_mat)
  !deallocate(Fp_mat)
  !deallocate(W)

  write(6,*)'end subroutine kmatrix'
End subroutine kmatrix


subroutine perturbation_99(max_open,K,K_perturbed,charge,energy,R0)
  use nrtype,only : i4b,dbl,dpc
  use solve, only : invert_real
  implicit none
  integer(kind=i4b),intent(in) :: max_open
  real(kind=dbl),intent(in) :: K(max_open,max_open),charge,energy,R0
  real(kind=dbl),intent(out) :: K_perturbed(max_open,max_open)
  real(kind=dbl),dimension(:,:),allocatable :: L,M,Gamma1,Gamma2,ss_c,sc_c,cs_c,cc_c,&
       & identity,K_pert_1,K_pert_2
  integer(kind=i4b) :: p,q


  allocate(ss_c(max_open,max_open))
  allocate(sc_c(max_open,max_open))
  allocate(cs_c(max_open,max_open))
  allocate(cc_c(max_open,max_open))
  allocate(Gamma1(max_open,max_open))
  allocate(Gamma2(max_open,max_open))
  allocate(identity(max_open,max_open))
  allocate(L(max_open,max_open))
  allocate(M(max_open,max_open))
  allocate(K_pert_1(max_open,max_open))
  allocate(K_pert_2(max_open,max_open))

  call Integral_coulomb(max_open,ss_c,sc_c,cs_c,cc_c)

  ! Perform the linear algebra
  identity=0.0d0
  do p=1,max_open
     do q=1,max_open
        if (p.eq.q) then
           identity(p,q)=1
        else
           identity(p,q)=0
        end if
     end do
  end do
  Gamma1=matmul(sc_c,K)
  Gamma2=matmul(cc_c,K)
  L(:,:)=-ss_c(:,:)+K(:,:)-Gamma1(:,:)
  M(:,:)=identity(:,:)+cs_c(:,:)+Gamma2(:,:)
  call invert_real(M,max_open,max_open)
  K_pert_1=matmul(L,M)
  Gamma1=matmul(K,cc_c)
  Gamma2=matmul(K,cs_c)
  L(:,:)=-ss_c(:,:)+K(:,:)-Gamma2(:,:)
  M(:,:)=identity(:,:)+sc_c(:,:)+Gamma1(:,:)
  call invert_real(M,max_open,max_open)
  K_pert_2=matmul(M,L)
  K_perturbed(:,:)=.5d0*(K_pert_1(:,:)+K_pert_2(:,:))


  deallocate(K_pert_1)
  deallocate(K_pert_2)
  deallocate(L)
  deallocate(M)
  deallocate(identity)
  deallocate(Gamma1)
  deallocate(Gamma2)
  deallocate(ss_c)
  deallocate(sc_c)
  deallocate(cs_c)
  deallocate(cc_c)

end subroutine perturbation_99


subroutine perturbation_87(max_open,R0,charge,energy,lmax,K_perturbed,Z_mat,b)
  use nrtype,only : i4b,dbl,dpc
  use solve, only : invert_real
  implicit none
  integer(kind=i4b),intent(in) :: max_open,lmax
  real(kind=dbl),intent(in) :: charge,energy,R0,Z_mat(max_open,max_open),b(max_open)
  real(kind=dbl),intent(out) :: K_perturbed(max_open,max_open)
  real(kind=dbl),dimension(:,:),allocatable :: T,M,Gamma1,Gamma2,ss_c,sc_c,cs_c,cc_c,&
       & identity,K_pert_1,K_pert_2,R
  real(kind=dbl),dimension(:),allocatable :: J,Jp,N,Np
  integer(kind=i4b) :: p,q,f

  allocate(ss_c(max_open,max_open))
  allocate(sc_c(max_open,max_open))
  allocate(cs_c(max_open,max_open))
  allocate(cc_c(max_open,max_open))
  allocate(Gamma1(max_open,max_open))
  allocate(Gamma2(max_open,max_open))
  allocate(identity(max_open,max_open))
  allocate(T(max_open,max_open))
  allocate(M(max_open,max_open))
  allocate(K_pert_1(max_open,max_open))
  allocate(K_pert_2(max_open,max_open))
  allocate(J(max_open))
  allocate(N(max_open))
  allocate(Jp(max_open))
  allocate(Np(max_open))
  allocate(R(max_open,max_open))

  write(26,*)max_open,charge,energy,R0
  write(27,*)J,Jp,N,Np
  write(28,*)ss_c,sc_c,cs_c,cc_c


  call Integral_Coulomb(max_open,ss_c,sc_c,cs_c,cc_c,charge,energy,R0,J,Jp,N,Np)

  R(:,:)=0.d0

  do p=1,max_open
     do q=1,max_open
        do f=1,max_open
           R(p,q)=R(p,q)+Z_mat(p,f)*Z_mat(q,f)/b(f)
        end do
     end do
  end do

  call write_simm(R,max_open,42)


  do p=1,max_open
     do q=1,max_open
        if (p.eq.q) then
           T(p,q) = N(p)
        else
           T(p,q) = 0.d0
        end if
     end do
  end do

  do p=1,max_open
     do q=1,max_open
        M(p,q) = T(p,q) - Np(p)*(R(p,q))
     end do
  end do

  call write_simm(M,max_open,43)


  write(6,*)'before invert_real'
  call invert_real(M,max_open,max_open)

  call write_simm(M,max_open,44)


  T=0.d0
  do p=1,max_open
     do q=1,max_open
        if (p.eq.q) then
           T(p,q) = J(p)
        else
           T(p,q) = 0.d0
        end if
     end do
  end do

  Gamma1=0.0d0

  do p=1,max_open
     do q=1,max_open
        Gamma1(p,q) = T(p,q) - Jp(p)*R(p,q)
     end do
  end do
  call write_simm(Gamma1,max_open,45)

  K_perturbed=matmul(Gamma1,M)

  call write_simm(K_perturbed,max_open,46)


  deallocate(J,Jp,N,Np)
  deallocate(K_pert_1)
  deallocate(K_pert_2)
  deallocate(T)
  deallocate(M)
  deallocate(identity)
  deallocate(Gamma1)
  deallocate(Gamma2)
  deallocate(ss_c)
  deallocate(sc_c)
  deallocate(cs_c)
  deallocate(cc_c)

end subroutine perturbation_87

subroutine get_qd(K,max_open,tan_qd,kappa)
  use nrtype,only : i4b,dbl,dpc,Pi
  implicit none
  integer(kind=i4b),intent(in) :: max_open
  real(kind=dbl), intent(in) :: K(max_open,max_open),kappa
  real(kind=dbl),intent(out) :: tan_qd(max_open)
  integer(kind=dbl) :: f

  call diag_symm(K,max_open,max_open,tan_qd)

  do f=1,max_open
     write(103,*)(kappa**2)/2,f,tan_qd(f)
  end do
  !Calculate the quantum defects
  write(6,*)'QD'
  do f=1,max_open
     tan_qd(f)=atan(tan_qd(f))/Pi
     write(476,*)(kappa**2)/2,f,tan_qd(f)
     write(6,*)(kappa**2)/2,f,tan_qd(f)
  end do

end subroutine get_qd

subroutine write_simm(K,max_open,number)
  use nrtype,only : i4b,dbl,dpc,Pi
  implicit none
  integer(kind=i4b),intent(in) :: max_open,number
  real(kind=dbl), intent(in) :: K(max_open,max_open)
  integer(kind=i4b) :: p,q

  do p=1,max_open
     do q=1,max_open
        write(number,*)K(p,q),p,q
     end do
  end do

end subroutine write_simm


