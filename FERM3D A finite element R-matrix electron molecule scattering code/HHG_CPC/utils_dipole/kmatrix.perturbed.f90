!!!*************************************************************
! 文件/File: kmatrix.perturbed.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: kmatrix.perturbed.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


subroutine kmatrix_perturbed(J,Jp,N,Np,charge,energy,H_vectors,KCC,Rfinal)
  use Solve,only : Z=>eigen_vector,b=>eigen_value,invert_real,invert_complex,diag_complex,diag_symm
  use Open_information,only : max_open,kappa,E,max_index,lmax,location,nE,S_eigenvec
  use Brep,only : R0,p_Patch,x_coord=>coord_x,y_coord=>coord_y,z_coord=>coord_z
  use nrtype, only :Pi,dbl,i4b,dpc
  use control , only : output_file,option_delay,option_wf,molecule
  use Matrices, only : RHS_co => d
  implicit none
  interface perturbation_87
     subroutine perturbation_87(max_open,R0,charge,energy,lmax,K_perturbed,Z_mat,b)
       use nrtype,only : i4b,dbl,dpc,Pi
       integer(kind=i4b),intent(in) :: max_open,lmax
       real(kind=dbl),intent(in) :: charge,energy,R0,Z_mat(max_open,max_open),b(max_open)
       real(kind=dbl),intent(out) :: K_perturbed(max_open,max_open)
     end subroutine perturbation_87
  end interface
  integer(kind=i4b) :: p,q,f,indice,alpha,enne,nClosed
  integer(kind=i4b),save :: counter=0
  real(kind=dbl),dimension(max_open) :: N2,J,Jp,N,Np,tan_qd
  real(kind=dbl),dimension(max_open) ::sigma_partial
  real(kind=dbl)::sigma_total,sigma_elastic,sigma_inelastic,sigma_elastic2,sigma_elastic3,sigma1,sigma2,Ene,eigenphase_sum,delta,eta,eigenphase_sum_2,charge,energy
  real(kind=dbl),dimension(:,:),allocatable :: R,K,M,T,mostro,astro,identity,I_mat,J_mat,F_mat,Fp_mat,Z_mat,cc,U_vec,K_perturbed
  real(kind=dbl) :: H_vectors(max_open,max_open),KCC(max_open,12)
  real(kind=dbl),dimension(:),allocatable :: W
  real(kind=dbl),parameter :: angstrom_atomicu_conversion=0.5291772083d0
  real(kind=dbl) :: x,y,zzz
  complex(kind=dpc),dimension(:,:),allocatable :: L,S,ident,probe,cc_complex,cc_complex2
  real(kind=dbl),dimension(:,:),allocatable,save ::coeff_projection
  complex(kind=dpc),dimension(:),allocatable :: diag_S
  complex(kind=dpc),parameter :: I=(0.0d0,1.0d0) 
  complex(kind=dpc) :: det,det2,argom,gamma_s,gammaln,cgamma
  complex(kind=dpc),parameter :: zI=(0.d0,1.d0)
  character*10 :: method='No perturbation'
  ! New variables x-form ---------------
  integer(kind=i4b) :: ll,mm,counter1,counter2
  real(kind=dbl) :: Rfinal
  complex(kind=dpc),dimension(:,:), allocatable :: T_mat 
  !-------------------------------------
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
  write(6,*)nClosed,p_Patch,max_index
  allocate(cc(nClosed,max_open))
  ! Debug
  !  cc= - matmul(RHS_co,Z)
  ! Call Blas 3 subroutine to perform this multiplication
  cc(:,:)=0.d0
  write(6,*)'dgemm0'
  if (option_wf.eq.'wf')call dgemm ( 'N', 'N', p_Patch-max_index, max_open, max_open, -1.d0, RHS_co, p_Patch-max_index, Z, max_open, &
       &                   0.d0, cc, p_Patch-max_index )

  cc=cc*sqrt(2.d0)

902 format(5e13.5)
  !**COMM**********
  N2=0.0d0
  R=0.0d0
  K=0.0d0
  M=0.0d0
  L=cmplx(0.0d0,0.0d0)
  S=cmplx(0.0d0,0.0d0)
  kappa=sqrt(2*E)

  do p=1,max_open
     do q=1,max_open
        do f=1,max_open
           N2(p)=N2(p)+Z(f,q)*Z(f,p)
           M(p,q)=M(p,q)+Z(f,q)*Z(f,p)
        end do
     end do
     if (N2(p).lt.0) then
        N2(p) = -N2(p)
     end if
  end do


  write(6,*)'Rmatrix'

  do f=1,max_open
     write(84,*)b(f)
     do q=1,max_open
        write(83,*)Z(f,q)
     end do
  end do

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
        write(65,*)R(p,q),p,q
     end do
  end do
  ! Invert (N-N'R)

  !------------------------------
  !X-form to spherical harmonic basis --> TEST
  if (counter==0) then
     allocate(coeff_projection(max_open,max_open))
     do p=1,max_open
        do q=1,max_open 
           read(30,*)coeff_projection(p,q)
        end do
     end do
     call invert_real(coeff_projection,max_open,max_open)
  end if
  counter=counter+1

  !------------------------
  ! MULTI
  ! Conversion from f form to u form of R-matrix:
  ! R(u)=ZR0*(b^-1*Z*R0+Z)^-1
  T(:,:)=0.d0
  do p=1,max_open
     do q=1,max_open
        T(p,q)=Z(p,q)*b(q)*R0+Z(p,q)
     end do
  end do
  M=T ! Z-prime in u-basis
  call invert_real(T,max_open,max_open)

  R=matmul(Z,T)
  R=R*R0 ! R-matrix in u-basis
  Z=Z*R0 ! Z x-form to u-form

  !------------------------------


  !--------------------------------
  ! Final transformation to short range R-matrix in real harmonics basis
  T=R
  T=matmul(T,transpose(coeff_projection/R0))  
  T=matmul(coeff_projection/R0,T)                   
  R=T
  !------------------------
  ! Test
  !write(300,*)'Rinitial'
  !do p=1,max_open
  !   write(300,300)(R(q,p),q=1,max_open)
  !end do

  !------------------------------

  Z=matmul(coeff_projection/R0,Z)  !MULTI: xform to real spherical harmonics
  M=matmul(coeff_projection/R0,M)  !MULTI: xform to real spherical harmonics

  if (molecule.eq.'neutral') then

  !------------------------
  ! Contraction of the R-matrix for dipole matching purposes to dipole harmonics
  T=0.d0
  T=matmul(R,transpose(H_vectors))
  R=matmul((H_vectors),T)
  M=matmul(H_vectors,M) ! x-form derivatives of eigenfunctions to dipole repr.
  Z=matmul(H_vectors,Z) ! x-form derivatives of eigenfunctions to dipole repr.
  !------------------------

  ! Propagation step (as in CW Clark PRA'79)
  allocate(F_mat(max_open,max_open))
  T(:,:)=0.d0
  mostro(:,:)=0.d0
  K=0.d0
  F_mat=0.d0

  do p=1,max_open
     T(p,p)=KCC(p,9)*KCC(p,8)-KCC(p,11)*KCC(p,6)
     mostro(p,p)=KCC(p,12)*KCC(p,5)-KCC(p,7)*KCC(p,10)
     K(p,p)=-KCC(p,9)*KCC(p,7)+KCC(p,11)*KCC(p,5)
     F_mat(p,p)=KCC(p,10)*KCC(p,8)-KCC(p,12)*KCC(p,6)
     !write(36,*)'T',T(p,p);   write(36,*)'mostro',mostro(p,p);   write(36,*)'K',K(p,p);   write(36,*)'F_mat',F_mat(p,p)
  end do
  !----------------------------------
  !T=matmul(T,Z)+matmul(K,M)
  !mostro=matmul(mostro,M)+matmul(F_mat,Z)  ! x-form to g,g' --> Works fine (slightly nonsymmetric)
  !----------------------------------
  T=matmul(T,R)
  F_mat=matmul(F_mat,R)
  T=matmul(T+K,M)
  mostro=matmul(mostro+F_mat,M) ! x-form to g',g' --> Works fine (symmetric)
  deallocate(F_mat)
  !----------------------------------

  call invert_real(mostro,max_open,max_open)
  R=matmul(T,mostro) ! THIS STILL DOES NOT WORK !!! WORKS WITHOUT DIPOLE MATCHING

  T=matmul(R,H_vectors)
  R=matmul(transpose(H_vectors),T)
  !------------------------------------------------
  ! Test
  !write(300,*)'Rfinal'
  !do p=1,max_open
  !   write(300,300)(R(q,p),q=1,max_open)
  !end do
300 format(100e15.8)
  end if ! End propagation for neutral molecules
  !------------------------------------------------



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

  write(6,*)'1st inversion'
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
  write(6,*)'before invert_real'
  call invert_real(M,max_open,max_open)
  call write_simm(M,max_open,53)

  write(6,*)'before mostro'
  mostro=matmul(T,M)

  ! K-matrix

  write(6,*)'invert_real'

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

  write(6,*)'K'

  K=matmul(Z,M)
  do p=1,max_open
     do q=1,max_open
        write(598,*)K(p,q),p,q
     end do
  end do
  !*COMM*****
  !*Diagonalize K-matrix to obtain the quantum defects
  write(6,*)'diag K'
  T=K


  call diag_symm(T,max_open,max_open,tan_qd)

  do f=1,max_open
     !write(103,*)(kappa**2)/2,f,tan_qd(f)
  end do
  !Calculate the quantum defects
  write(6,*)'QD'
  do f=1,max_open
     tan_qd(f)=atan(tan_qd(f))/Pi
     write(473,*)(kappa**2)/2,f,tan_qd(f)
     write(6,*)(kappa**2)/2,f,tan_qd(f)
  end do

  !-------------------------------------------------------------------------------------
  ! Perturbative (multipole) part
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
  !write(472,*)max_open
  do p=1,max_open
     do q=1,max_open
        !write(472,*)T(p,q)
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
  call invert_real(I_mat,max_open,max_open)


  ! Not good when dipole matching, can't do contraction on K-matrix
  !  K=matmul(J_mat,I_mat)

  !------------------------
  ! Contraction of the K-matrix for dipole matching purposes (does not work on
  ! K-matrix)
  !T=0.d0
  !T=matmul(K,H_vectors)
  !K=matmul(transpose(H_vectors),T)
  !------------------------

  do p=1,max_open
     do q=1,max_open
        write(600,*)K(p,q),p,q
     end do
  end do
  T=K
  call diag_symm(T,max_open,max_open,tan_qd)



  sigma_elastic=0.d0
  do f=1,max_open
     tan_qd(f)=atan(tan_qd(f))/Pi
     write(490,*)(kappa**2)/2,f,tan_qd(f)
     sigma_elastic=sigma_elastic+Pi/(kappa*kappa)*4*sin(Pi*tan_qd(f))**2
  end do
  !write(701,*)sigma_elastic,kappa**2/2

  !S_eigenvec(:,:,nE)=cmplx(T(:,:),0.d0) ! If need to put K-eigenvectors on S-solution
  if (option_wf.eq.'wf') then
     if (.not.option_delay) then
        !*******************************************************
        ! ** Transformations to plot the wavefunction -- According to Chris'
        ! ** suggestions for ions
        ! ** Transformation to collisional eigenchannels

        !  do p=1,max_open
        !     do q=1,max_open
        !        T(p,q)=T(p,q)*cos(Pi*tan_qd(q))
        !     end do
        !  end do

        !  cc=matmul(cc,I_mat)
        !  cc=matmul(cc,T)
        ! Call Blas 3 subroutine to perform this multiplication
        cc(:,:)=0.d0
        !call dgemm ( 'N', 'N', p_Patch-max_index, max_open, max_open, 1.d0, cc, p_Patch-max_index, I_mat, max_open, &
        !     &                   0.d0, cc, p_Patch-max_index ) ! Wrong at this point, need to correct
        write(6,*)'dgemm1',p_Patch-max_index, max_open
        call dgemm ( 'N', 'N', p_Patch-max_index, max_open, max_open, 1.d0, RHS_co, p_Patch-max_index, T, max_open, &
             &                   0.d0, cc, p_Patch-max_index )
        !  do f=1,max_open
        !     q=0
        !     do p=1,nClosed-max_index,8
        !        q=q+1
        !        if (f.eq.location) write(10001,902)x_coord(q),y_coord(q),z_coord(q),cc(p,f)
        !     end do
        !  end do
903     format(4e13.5,I)
     else  ! if no K-matrix matching in time-delay calc

        !-----------------------------------------
        ! Transformation to S matrix eigenchannels, plot the density
        L(:,:)=(M(:,:)+cmplx(0.d0,1.d0)*J_mat(:,:))
        call invert_complex(L,max_open,max_open)
        ! L is I_mat+I*J_mat
        !  cc=matmul(cc,L)
        !  cc=matmul(cc,T)
        ! Call Blas 3 subroutine to perform this multiplication
        cc(:,:)=0.d0
        allocate(cc_complex((p_Patch-max_index),max_open))
        allocate(cc_complex2((p_Patch-max_index),max_open))

        allocate(probe(max_open,max_open))
        probe(:,:)=S_eigenvec(1:max_open,1:max_open,nE)
        cc_complex(:,:)=cmplx(0.d0,0.d0)

        !  do p=1,max_open
        !     do q=1,max_open
        !        probe(p,q)=probe(p,q)*cos(Pi*tan_qd(q))
        !     end do
        !  end do
        cc_complex2(:,:)=cmplx(RHS_co(:,:),0.d0)
        write(6,*)'zgemm1',p_Patch-max_index, max_open,p_Patch,max_index
        call zgemm ( 'N', 'N', p_Patch-max_index, max_open, max_open, (1.d0,0.d0), cc_complex2, p_Patch-max_index, L, max_open, &
             &                   (0.d0,0.d0), cc_complex, p_Patch-max_index )
        !cc_complex2(:,:)=cmplx(RHS_co(:,:),0.d0)
        cc_complex2(:,:)=cc_complex(:,:)
        cc_complex(:,:)=cmplx(0.d0,0.d0)
        write(6,*)'zgemm2'
        call zgemm ( 'N', 'N', p_Patch-max_index, max_open, max_open, (1.d0,0.d0), cc_complex2, p_Patch-max_index, probe, max_open, &
             &                   (0.d0,0.d0), cc_complex, p_Patch-max_index )
        write(6,*)'after zgemm'
        ! Rescaling for phase to get dominant real part
        indice=maxloc(abs(cc_complex(:,location)),1)
        det=0.d0
        if (abs(dble(cc_complex(indice,location))).gt.1.d-12) det=aimag(cc_complex(indice,location)) / dble(cc_complex(indice,location))
        cc_complex(:,location)=cc_complex(:,location)*exp(-zI*atan( dble(det) ))
        !----------------------------------------------
        do f=1,max_open
           q=0
           do p=1,nClosed-max_index,8
              q=q+1
              if (f.eq.location)    then
                 write(11001,902)x_coord(q),y_coord(q),z_coord(q),dble(cc_complex(p,f))
                 write(12001,902)x_coord(q),y_coord(q),z_coord(q),aimag(cc_complex(p,f))
                 write(13001,902)x_coord(q),y_coord(q),z_coord(q),abs(cc_complex(p,f))
              end if
              write(10000+f,902)x_coord(q),y_coord(q),z_coord(q),dble(cc_complex(p,f)),aimag(cc_complex(p,f))
           end do
        end do
        deallocate(cc_complex)
        deallocate(cc_complex2)

        deallocate(probe)
     end if
  end if
  !-----------------------------------------
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

  eigenphase_sum=0.d0
  eigenphase_sum_2=0.d0
  do f=1,max_open
     eigenphase_sum=eigenphase_sum+tan_qd(f)
     eigenphase_sum_2=eigenphase_sum_2+abs(tan_qd(f))
  end do
  write(475,*)eigenphase_sum,(kappa**2)/2
  write(469,*)eigenphase_sum_2*Pi,(kappa**2)/2

  !*COMM*******
  !*Calculation of the energies

  !do alpha=1,max_open
  !   do enne=1,max_open
  !      Ene=-1.d0/(2.d0*(enne-tan_qd(alpha))**2)
  !      write(473,*)Ene,enne,alpha
  !   end do
  !end do

  !*COMM*****

  write(6,*)'2nd inversion'
  L=cmplx(0.0d0,0.0d0)

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

  allocate(probe(max_open,max_open))
  allocate(ident(max_open,max_open))
  probe=cmplx(0.0d0,0.0d0)
  probe=L
  ident=cmplx(0.0d0,0.0d0)

  call invert_complex(L,max_open,max_open)
  write(6,*)'after invert complex'
  ident=matmul(probe,L)
  write(6,*)'after matmul'

  S=cmplx(0.0d0,0.0d0)
  S=matmul(cmplx(identity,K),L)

  !------------------------------
  !X-form to spherical harmonic basis --> TEST
  !allocate(coeff_projection(max_open,max_open))
  !  do p=1,max_open
  !    do q=1,max_open 
  !      read(30,*)coeff_projection(p,q)
  !    end do
  !  end do
  !----------------- X-form to complex spherical harmonics --------------
!!$ Transform to spherical harmonics
  ! Numerical_H=T^-1 * C^-1 * Ylm
  ! Construct T matrix
  allocate(T_mat(max_open,max_open))
  T_mat(:,:)=(0.d0,0.d0)
  do ll=0,lmax-1
     do mm=-ll,ll
        counter1=ll**2+1+(ll+abs(mm))
        counter2=ll**2+1+(ll-abs(mm))
        if (mm.gt.0) then
           T_mat(counter1,counter2)=1.d0/sqrt(2.d0) ! +mm
           T_mat(counter1,counter1)=1.d0/sqrt(2.d0) ! -mm

        else if (mm.lt.0) then
           T_mat(counter2,counter1)=1.d0/sqrt(2.d0)/zI ! +mm
           T_mat(counter2,counter2)=-1.d0/sqrt(2.d0)/zI ! -mm

        else if (mm.eq.0) then
           T_mat(counter1,counter2)=1.d0
        end if
        !write(6,*)counter1,counter2,T_mat(counter1,counter2)

     end do
  end do
  !T_mat(:,:)=(0.d0,0.d0) ! Use as a test, to see only coeff_projection
  do ll=1,max_open
     do mm=1,max_open
        !if (ll.eq.mm) T_mat(ll,mm)=(1.d0,0.d0)
        write(32,*)T_mat(ll,mm),real(coeff_projection(ll,mm))
     end do
  end do

  ! T_mat(:,:)=matmul(coeff_projection(:,:),T_mat(:,:))

  call invert_complex(T_mat,max_open,max_open)
  write(32,*)T_mat(:,:)
  ! coeff_projection(:,:)=T_mat(:,:) ! Pass the final transformation to the main routine

  !---------------- End modified x-form -----------------------
  write(6,*)'after invert complex xform' 
  write(6,*)'counter=',counter
  ident=S
  ident=matmul(ident,transpose(conjg(T_mat)))
  ident=matmul(T_mat,ident)
  S=ident
  deallocate(T_mat)

  !deallocate(coeff_projection)
  !------------------------------
  do p=1,max_open
     do q=1,max_open 
        write(538,*)S(p,q),p,q
     end do
  end do

  ! Diagonalization of the S matrix
  write(6,*)'before alloc'
  allocate(diag_S(max_open))
  write(6,*)'after alloc'
  !***COMM***
  !*Diagonalize S-matrix

  call diag_complex(S,L,max_open,max_open,diag_S)
  write(6,*)'after diag_complex'
  !  do p=1, max_open
  !     do q=1,max_open
  !        write(399,*)S(p,q),p,q
  !     end do
  !    write(398,*)diag_S(p),p
  !  end do
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
  write(6,*)'eigenphase sum'
  eigenphase_sum=0.d0
  do p=1,max_open
     !     delta=asin(.25d0*sqrt(real((diag_S(p)-1)*conjg(diag_S(p)-1))))
     delta=.5d0*acos(1.-real((diag_S(p)-1)*conjg(diag_S(p)-1))/2.d0)
     eigenphase_sum=eigenphase_sum+delta
  end do
  write(478,*)eigenphase_sum,(kappa**2)/2
  !  eigenphase_sum=0.d0
  !  do p=1,max_open
  !     if (abs(diag_S(p)).ge.1.d-5) then
  !        delta=dble(-.5d0*I*(log(diag_S(p))))
  !write(402,*)delta,p
  !        eigenphase_sum=eigenphase_sum+abs(delta)
  !     end if
  !  end do
  !  write(384,*)eigenphase_sum,(kappa**2)/2
  !  det=(1.d0,0.d0)

  !  do p=1,max_open
  !     det=det*S(p,p)
  !  end do
  !  delta=asin(.25d0*sqrt(real((det-1)*conjg(det-1))))


  !  write(479,*)delta,(kappa**2)/2

  !  sigma_total=0.0d0
  sigma_elastic=0.0d0
  !  sigma_elastic2=0.0d0
  !  sigma_elastic3=0.0d0
  !  sigma_inelastic=0.0d0
  write(6,*)'cross sections'
  !deallocate(S)
  ! Cross sections
  ! Total
  !  do p=1,max_open
  !     sigma_total=sigma_total+  Pi / (kappa*kappa) * (1- real(S(p,p)))
  !  end do
  !  write(380,*)sigma_total, kappa**2 / 2
  ! Elastic
  do p=1,max_open
     !     sigma_elastic=sigma_elastic+Pi/(kappa*kappa)*(S(p,p)-1)*conjg(S(p,p)-1) 
     sigma_elastic=sigma_elastic+Pi/(kappa*kappa)*(diag_S(p)-1)*conjg(diag_S(p)-1) 

     !     sigma_elastic2=sigma_elastic2+Pi/(kappa*kappa)*(S(p,p)-1)*conjg(S(p,p)-1)
  end do
  !  write(983,*)sigma_elastic, kappa**2 / 2
  !write(output_file,*)sigma_elastic, kappa**2 / 2
  !Cross section in square angstroms
  !  write(985,*)sigma_elastic*angstrom_atomicu_conversion**2,kappa**2/2
  write(6,*)'cross section=',sigma_elastic, kappa**2 / 2.d0
  write(1280,*)kappa**2/2.d0,sigma_elastic

!!!$$$ Rewrite partial cross sections
  !do p=0,lmax-1
  !  sigma_partial(p)=(2*p+1)/(4*kappa*kappa)*(S(p,p)-1)*conjg(S(p,p)-1)
  !write(400+p,*)sigma_partial(p),kappa**2 /2
  !write(600+p,*)sigma_partial(p),kappa
  !end do
!!!$$$
  write(6,*)'before deallocations'
  deallocate(cc)
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

  write(6,*)'after deallocations'
End subroutine kmatrix_perturbed


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
