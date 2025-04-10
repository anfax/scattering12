!!!*************************************************************
! 文件/File: R_matrix_post_production.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: R_matrix_post_production.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


module alpha
contains
  subroutine R_matrix_post_production
    use nrtype
    use control
    use Brep, only : R0,coord_x,coord_y,coord_z,p_Patch
    use Solve
    use Open_information,only : kappa,E,max_open,lmax,max_index,location,nE,S_eigenvec,max_l
    use Matrices, only : RHS_co => d
    use Hydrogen

    implicit none
    !Vars for channel_computation
    integer(kind=i4b) :: nEmax,lds,i,j,k,nClosed,ios_DCS
    integer(kind=i4b),dimension(:),allocatable :: max_open_array,l_max
    integer(kind=i4b),dimension(:),allocatable :: loc
    real(kind=dbl) :: Emax,Emin,Delta_E
    complex(kind=dpc),dimension(:,:,:),allocatable :: S
    complex(kind=dpc),dimension(:,:),allocatable :: S_current,S1,S2
    complex(kind=dpc),dimension(:),allocatable :: eigen,eigen1,eigen2
    real(kind=dbl), dimension(:),allocatable :: b,H_values
    real(kind=dbl), dimension(:,:),allocatable :: qd,H_vectors,KCC
    real(kind=dbl), dimension(:,:,:),allocatable :: U
    real(kind=dbl) :: absolute,temp,temp2,dipole(5)=0.d0,dipole_temp(5)=0.d0
    integer(kind=i4b), parameter :: nx=10,nw=12
    character(10) :: option_DCS='no'
    character(15) :: tempchar='temp'
    logical :: rotation
    !Input all data

    ! Read "input_grid" and "input_control"

    !--------------------
    ! Open files
    open(unit=nx,file='input_grid.dat',status='old',action='read')
    open(unit=nw,file='input_control.dat',status='old',action='read')
    !--------------------
    ! control
    read(nw,*)output_file
    read(nw,*)tempchar,option_wf
    read(nw,*)tempchar,max_l
    read(nw,*)molecule,Z
    write(6,*)Z
    read(nw,*)
    read(nw,*)
    read(nw,*)
    read(nw,*)
    read(nw,*)
    read(nw,*)Emax
    read(nw,*)Emin
    read(nw,*)nEmax
    read(nw,*)
    read(nw,*)rotation
    read(nw,*)
    read(nw,*)
    read(nw,*)
    read(nw,*,IOSTAT=ios_DCS)option_DCS,dipole(1),dipole(2),dipole(3),dipole(4)
 !write to std output the input flags
  write(6,*)'-----------------------------------------'
  write(6,*)'CALCULATION FLAGS'
  write(6,*)'output file= ',output_file
  write(6,*)'calculation type= ',tempchar,'wf= ',option_wf
  write(6,*)'molecule type= ',molecule,'charge= ',Z
  write(6,*)'rotation= ',rotation
  write(6,*)'energy= ',Emin,Emax,nEmax

  write(6,*)'-----------------------------------------'


    ! Take into account rotation of the molecule
    if (rotation) then
         dipole(5)=-1.d0
    else
         dipole(5)=1.d0
    end if 

    write(6,*)'options',option_wf,option_DCS,dipole
    option_delay=.TRUE.
9000 format(2a15)
9001 format(1a15,1f13.5)
    !--------------------
    !--------------------
    ! grid
    read(nx,*)R0

    !--------------------
    ! Close files
    close(nx)
    close(nw)
    !--------------------

    allocate(max_open_array(nEmax),l_max(nEmax),loc(nEmax))
    call channel_calculation(Emax,Emin,nEmax,max_open_array,l_max,Delta_E)
    ! Calc max number of channels at each energy
    max_open=max_open_array(nEmax)
    !------------------------------------------------
    ! Dipole part
    !allocate(H_values(max_open))
    !allocate(H_vectors(max_open,max_open))
    !dipole=4.7d0
    !call Hydrogen_perturbed(dipole,H_vectors,H_values,max_open)
    !stop
    !------------------------------------------------
    allocate(S(max_open,max_open,nEmax))
    allocate(S_eigenvec(max_open,max_open,nEmax))
    allocate(U(max_open,max_open,nEmax))
    allocate(qd(max_open,nEmax))
    S(:,:,:)=cmplx(0.d0,0.d0)
    U(:,:,:)=0.d0
    qd(:,:)=0.d0
    write(6,*)nEmax*max_open*max_open,max_open_array(nEmax),nEmax
    do nE=1,nEmax
       ! Read the S(537)
       max_open=max_open_array(nE)
       write(6,*)'nE',max_open,nE
       read(472,*)
       do j=1,max_open
          do k=1,max_open
             if (option_delay) then
                read(536,*)S(j,k,nE)
             else
                read(536,*)S(j,k,nE)
                read(472,*)U(j,k,nE)
             end if
          end do
          read(474,*)temp,qd(j,nE)
       end do
    end do
    write(6,*)'end reading',nEmax
    max_open=max_open_array(nEmax)
    if (option_delay) then
       call time_delay(nEmax,max_open, max_open,Emax,Emin,loc,max_open_array,l_max,S)
    else
       call time_delay2(nEmax,max_open, max_open,Emax,Emin,loc,max_open_array,l_max,U,qd)
    end if
    ! Differential cross sections
    if (option_DCS.eq.'yes') then
      call differential(nEmax,max_open, max_open,Emax,Emin,loc,max_open_array,l_max,U,qd,S)
    end if
    do nE=1,nEmax
       max_open=max_open_array(nE)
       write(6,*)'max_open=',max_open,Emin,Emax,Delta_E
       lmax=l_max(nE)
       E=Emin+nE*Delta_E
       kappa=sqrt(2.d0*E)
       ! Read the RHC_co(?), Z(?)-->change name it's called like the charge, b(?),[x_coord,y_coord,z_coord](8)
       read(66,*)nClosed,p_Patch,max_index
       write(6,*)'max_index=',max_index,nClosed,E
       allocate(RHS_co(nClosed,max_open))
       allocate(coord_x(nClosed))
       allocate(coord_y(nClosed))
       allocate(coord_z(nClosed))
       allocate(eigen_value(max_open))
       allocate(eigen_vector(max_open,max_open))
       write(6,*)'max_index=',max_index,nClosed
       location=loc(nE)
       call input_matrices(nClosed,max_open,max_index,p_Patch,RHS_co,eigen_vector,eigen_value,coord_x,coord_y,coord_z)
write(6,*)'dipole=',dipole,Z
if ((molecule.eq.'neutral')) then
    !------------------------------------------------
    ! Dipole part
    allocate(H_values(max_open))
    allocate(H_vectors(max_open,max_open))
    dipole_temp=dipole
    call Hydrogen_perturbed(dipole_temp(1:5),H_vectors,H_values,max_open)
    !------------------------------------------------

allocate(KCC(max_open,12))
call continuum(Z,R0,H_vectors,max_open,H_values,KCC,E)
deallocate(KCC)
else
    allocate(H_values(max_open))
    allocate(H_vectors(max_open,max_open))
H_vectors=0.d0
H_values=0.d0
do i=1,max_open
H_vectors(i,i)=1.d0
end do
       call cross_section(H_vectors)
end if
       deallocate(RHS_co)
       deallocate(coord_x,coord_y,coord_z)
       deallocate(eigen_value,eigen_vector)
    deallocate(H_vectors,H_values)

    end do


!----------------------------
!    allocate(S(max_open,max_open,nEmax))
!    allocate(S_eigenvec(max_open,max_open,nEmax))
!    allocate(U(max_open,max_open,nEmax))
!    allocate(qd(max_open,nEmax))
!    S(:,:,:)=cmplx(0.d0,0.d0)
!    U(:,:,:)=0.d0
!    qd(:,:)=0.d0
!    write(6,*)nEmax*max_open*max_open,max_open_array(nEmax),nEmax
!    rewind(538)
!    do nE=1,nEmax
       ! Read the S(537)
!       max_open=max_open_array(nE)
!       write(6,*)'nE',max_open,nE
!       read(472,*)
!       do j=1,max_open
!          do k=1,max_open
!             if (option_delay) then
!                read(538,*)S(j,k,nE)
!             else
!                read(538,*)S(j,k,nE)
!                read(472,*)U(j,k,nE)
!             end if
!          end do
!          read(474,*)temp,temp2,qd(j,nE)
!       end do
!    end do
!    write(6,*)'end reading',nEmax
!    max_open=max_open_array(nEmax)
!    if (option_delay) then
!       call time_delay(nEmax,max_open, max_open,Emax,Emin,loc,max_open_array,l_max,S)
!    else
!       call time_delay2(nEmax,max_open, max_open,Emax,Emin,loc,max_open_array,l_max,U,qd)
!    end if
    ! Differential cross sections
!    if (option_DCS.eq.'yes') then
!      call differential(nEmax,max_open, max_open,Emax,Emin,loc,max_open_array,l_max,U,qd,S)
!    end if


!----------------------------






    deallocate(U,qd)
    deallocate(S)
    deallocate(loc)
    deallocate(S_eigenvec)

  end subroutine R_matrix_post_production

  subroutine channel_calculation(Emax,Emin,nEmax,max_open_array,l_max,Delta_E)
    use nrtype,only : i4b, dbl, dpc
    use Brep, only : R0
    use Open_information, only : max_l 
    implicit none
    integer(kind=i4b) :: nEmax
    real(kind=dbl) :: Emax,Emin,Delta_E,nE,max_open,li
    integer(kind=i4b),dimension(nEmax) :: max_open_array,l_max
    logical,dimension(:,:),allocatable :: opench
    real(kind=dbl), dimension(:,:),allocatable :: Ekin,T
    real(kind=dbl), dimension(:), allocatable :: Energy
    real(kind=dbl) :: E
    integer(kind=i4b),parameter :: mmax=100
    Delta_E=(Emax-Emin)/dble(nEmax)
    nEmax=anint((Emax-Emin)/Delta_E)
    !----------------------------------------------------------------------------
    !Allocation
    allocate(opench(mmax+1,nEmax))
    allocate(Ekin(mmax+1,nEmax))
    allocate(T(mmax+1,nEmax))
    allocate(Energy(nEmax))
    !----------------------------------------------------------------------------
    !----------------------------------------------------------------------------
    ! WARNING
    ! The open channels should not exceed the number of the open functions
    !-------------
    !Initialization
    opench(:,:)=.FALSE.
    Ekin(:,:)=0.d0
    T(:,:)=0.d0
    Energy(:)=0.d0
    max_open_array(:)=0
    l_max(:)=0
    !-------------

    do nE=1,nEmax
       E=Emin+nE*Delta_E
       Energy(nE)=E
       max_open=0
       do li=0,mmax
          T(li+1,nE)=li*(li+1)/R0**2.0d0
          Ekin(li+1,nE)=E-T(li+1,nE)
          if (Ekin(li+1,nE) .gt. 0.0d0) then
             opench(li+1,nE)=.TRUE.
             max_open=max_open+1
          else
             opench(li+1,nE)=.FALSE.
          end if
       end do
!       l_max(nE)=max_open+3
!       max_open_array(nE)=l_max(nE)*l_max(nE)
     l_max(nE)=max_l
     max_open_array(nE)=l_max(nE)**2
       !write(6,*)'l_max',myrank,l_max(nE),max_open_array(nE)
    end do
    deallocate(opench)
    deallocate(Ekin)
    deallocate(T)
    deallocate(Energy)
  end subroutine channel_calculation

  subroutine time_delay(nEmax,max_open,lds,Emax,Emin,loc,max_open_array,l_max,S)
    use nrtype, only : dpc,dbl,i4b,Pi
    use Open_information, only : S_eigenvec
    use Solve
    implicit none
    integer(kind=i4b) :: nEmax,max_open,lds,nE,lmax,max_open_top,i,max_l,l,m,channel
    integer(kind=i4b),dimension(nEmax) :: loc
    real(kind=dbl) :: Emax,Emin,Delta_E,E,eigenphase_sum,result
    integer(kind=i4b),dimension(nEmax) :: max_open_array,l_max
    complex(kind=dpc),dimension(max_open,max_open,nEmax) :: S
    complex(kind=dpc),dimension(:,:),allocatable :: S_current,S1,S2,Identity,dummy_mat
    complex(kind=dpc),dimension(:),allocatable :: eigen,eigen1,eigen2
    complex(kind=dpc), parameter :: zI=(0.d0,1.d0)
    integer(kind=i4b) :: j
    real(kind=dbl),dimension(:),allocatable :: absolute
    Delta_E=(Emax-Emin)/dble(nEmax)
    write(6,*)'inside time_del',nEmax,max_open
    write(6,*)max_open_array(:),l_max(:)
    max_open_top=max_open_array(nEmax)
    loc(1)=1
    loc(nEmax)=1
    do nE=2,nEmax-1
       max_open=max_open_array(nE)
       lmax=l_max(nE)
       E=Emin+nE*Delta_E
       write(6,*)max_open,lmax,E,Emin,Emax,Delta_E
       allocate(Identity(max_open_top,max_open_top))
       allocate(eigen(max_open_top),eigen1(max_open_top),eigen2(max_open_top),absolute(max_open_top))
       allocate(S_current(max_open_top,max_open_top))
       allocate(S1(max_open_top,max_open_top))
       allocate(S2(max_open_top,max_open_top))
       allocate(dummy_mat(max_open_top,max_open_top))

       Identity(:,:)=(0.d0,0.d0)
       do j=1,max_open_top
          Identity(j,j)=(1.d0,0.d0)
       end do
       S_current(:,:)=cmplx(0.d0,0.d0)
       S1(:,:)=cmplx(0.d0,0.d0)
       S2(:,:)=cmplx(0.d0,0.d0)

       S_current(:,:)=S(:,:,nE)
       S1(:,:)=S(:,:,nE-1)
       S2(:,:)=S(:,:,nE+1)

       do i=max_open+1, max_open_top
          S_current(i,i)=cmplx(1.d0,0.d0)
       end do
       do i=max_open_array(nE-1)+1, max_open_top
          S1(i,i)=cmplx(1.d0,0.d0)
       end do
       do i=max_open_array(nE+1)+1, max_open_top
          S2(i,i)=cmplx(1.d0,0.d0)
       end do

       write(6,*)'after allocs'
       !------------------------------
       S2(:,:)= (S2(:,:)-S1(:,:))/(2.d0*Delta_E)
       S2=zI*matmul(S_current,conjg(transpose(S2)))

       S1(:,:)=matmul(S2,transpose(conjg(S2)))
! Check S1
!write(34,*)'S_1=',S1(:,:)
! Make S2 (Q) hermitian
       S1(:,:)=(S2(:,:)+transpose(conjg(S2(:,:))))/2.d0
       S2=S1

       do i=1,max_open_top
          write(43,*)abs(S2(i,i)),nE,i
       end do
       call diag_complex(S2,Identity,max_open_top,max_open_top,eigen2)
       do i=1,max_open_top
          write(44,*)abs(eigen2(i)),nE,i
          write(42,*)eigen2(i),nE,i
       end do
       loc(nE)=maxloc(abs(eigen2(:)),1) ! ARE YOU SURE IT WORKS?
       write(43,*)'vector=',loc(nE),S2(:,loc(nE))
       Identity(:,:)=S2(:,:)


       write(43,*)'vector_real=',loc(nE),S2(:,loc(nE))

       S_eigenvec(:,:,nE)=S2(:,:)
       do i=max_open_array(nE+1), max_open_top
          S_eigenvec(i,i,nE)=cmplx(0.d0,0.d0)
       end do

       !-------------------------------
       deallocate(Identity)

       write(45,*)loc(nE),E
       ! Get the eigenphase sum, get resonance parameters
       call diag_complex(S_current,Identity,max_open_top,max_open_top,eigen)
       eigen(:)=log((eigen(:)))/(2.d0*zI)

       eigenphase_sum=0.d0
       eigenphase_sum=sum(real(eigen(:)))
       write(46,*)E,eigenphase_sum/Pi
       write(34,*)E,sum(abs(eigen2(:)))
       !--------------------------------------------------
       !--------------------------------------------------
       ! Get dominant partial wave contribution
       max_l=maxloc(real(S_eigenvec(:,loc(nE),nE)*S_eigenvec(:,loc(nE),nE)),1)
       write(41,*)max_l,nE,S_eigenvec(:,loc(nE),nE)*S_eigenvec(:,loc(nE),nE)
       ! sum the l components
       do l=0,lmax
        result =0.d0
        do m=1, 2*l+1
        channel=(l)**2+m
           result= result+real(S_eigenvec(channel,loc(nE),nE)*conjg(S_eigenvec(channel,loc(nE),nE)))
        write(51,*)maxloc(real(S_eigenvec(:,channel,nE)*S_eigenvec(:,channel,nE)),1),nE,real(eigen2(channel))
        end do
       write(50,*)result, l, nE
       end do
       !--------------------------------------------------

       deallocate(eigen,eigen1,eigen2,absolute)
       deallocate(S_current)
       deallocate(S1)
       deallocate(S2)
       deallocate(dummy_mat)
    end do

  end subroutine time_delay

  subroutine  time_delay2(nEmax,max_open,lds,Emax,Emin,loc,max_open_array,l_max,U,qd)
    use nrtype, only : dpc,dbl,i4b,Pi
    use Solve
    implicit none
    integer(kind=i4b) :: nEmax,max_open,lds,nE,lmax,max_open_top,i
    integer(kind=i4b),dimension(nEmax) :: loc
    real(kind=dbl) :: Emax,Emin,Delta_E,E,eigenphase_sum,max_l,result
    integer(kind=i4b),dimension(nEmax) :: max_open_array,l_max
    real(kind=dbl),dimension(max_open,max_open,nEmax) :: U
    real(kind=dbl),dimension(max_open,max_open,nEmax) :: U_copy
    real(kind=dbl),dimension(max_open,nEmax) :: qd
    real(kind=dbl),dimension(:),allocatable :: dot1,dot2
    integer(kind=i4b),dimension(:),allocatable :: top1,top2
    !  complex(kind=dpc),dimension(max_open,max_open,nEmax) :: S
    complex(kind=dpc),dimension(:,:),allocatable :: S_current,S1,S2,Identity
    complex(kind=dpc),dimension(:),allocatable :: eigen,eigen1,eigen2
    complex(kind=dpc),dimension(:),allocatable :: eigen_copy1,eigen_copy2
    complex(kind=dpc), parameter :: zI=(0.d0,1.d0)
    integer(kind=i4b) :: j,channel,l,m
    real(kind=dbl),dimension(:),allocatable :: absolute

U_copy(:,:,:)=0.d0
    Delta_E=(Emax-Emin)/dble(nEmax)
    write(6,*)'inside time_del2',nEmax,max_open
    write(6,*)max_open_array(:),l_max(:)
    max_open_top=max_open_array(nEmax)
    loc(1)=1
    loc(nEmax)=1
    do nE=2,nEmax-1
       max_open=max_open_array(nE)
       lmax=l_max(nE)
       E=Emin+nE*Delta_E
       !write(6,*)max_open,lmax,E,Emin,Emax,Delta_E
       allocate(Identity(max_open_top,max_open_top))
       allocate(eigen(max_open_top),eigen1(max_open_top),eigen2(max_open_top),absolute(max_open_top))
       allocate(S_current(max_open_top,max_open_top))
       allocate(S1(max_open_top,max_open_top))
       allocate(S2(max_open_top,max_open_top))

       Identity(:,:)=(0.d0,0.d0)
       do j=1,max_open_top
          Identity(j,j)=(1.d0,0.d0)
       end do
       !S_current(:,:)=S(:,:,nE)
       !S1(:,:)=S(:,:,nE-1)
       !S2(:,:)=S(:,:,nE+1)
       eigen(:)=exp(2.d0*zI*Pi*qd(:,nE))
       eigen1(:)=exp(2.d0*zI*Pi*qd(:,nE-1))
       eigen2(:)=exp(2.d0*zI*Pi*qd(:,nE+1))

       do i=max_open, max_open_top
          !S_current(i,i)=cmplx(1.d0,0.d0)
          eigen(i)=cmplx(1.d0,0.d0)
       end do
       do i=max_open_array(nE-1)+1, max_open_top
          !S1(i,i)=cmplx(1.d0,0.d0)
          eigen1(i)=cmplx(1.d0,0.d0)
       end do
       do i=max_open_array(nE+1)+1, max_open_top
          !S2(i,i)=cmplx(1.d0,0.d0)
          eigen2(i)=cmplx(1.d0,0.d0)
       end do
       write(6,*)'after allocs'
       !--------------------------------------
       ! Reordering
       allocate(eigen_copy1(max_open_top))
       allocate(eigen_copy2(max_open_top))
       allocate(top1(max_open_top))
       allocate(top2(max_open_top))
       allocate(dot1(max_open_top))
       allocate(dot2(max_open_top))

       eigen_copy1(:)=eigen1(:)
       eigen_copy2(:)=eigen2(:)
       do j=1,max_open_top
          do i=1,max_open_top
             dot1(i)=dot_product(U(:,i,nE-1),U(:,j,nE))
             dot2(i)=dot_product(U(:,i,nE+1),U(:,j,nE))
             !write(36,*)dot1(i),dot2(i),i,j
          end do
          top1(j)=maxloc(abs(dot1(:)),1)
          top2(j)=maxloc(abs(dot2(:)),1)
          write(37,*)j,top1(j),top2(j),E
       end do
       do i=1,max_open_top
          eigen1(i)=eigen_copy1(top1(i))
          eigen2(i)=eigen_copy2(top2(i))
          U_copy(:,i,nE-1)=U(:,top1(i),nE-1)
          U_copy(:,i,nE+1)=U(:,top2(i),nE+1)

          if (i.gt.max_open) then
             eigen1(i)=cmplx(1.d0,0.d0)
             eigen2(i)=cmplx(1.d0,0.d0)
          end if
       end do
       !--------------------------------------
       !------------------------------
       ! 1

       do i=1,max_open
          write(48,*)(eigen1(i)),(eigen2(i)),i,E,qd(top1(i),nE-1),qd(top2(i),nE+1)
       end do
U_copy(:,:,nE)=(U_copy(:,:,nE+1)+U_copy(:,:,nE-1))/(2.d0*Delta_E)*U(:,:,nE)
write(40,*)nE
write(40,*)U_copy(:,:,nE)
       eigen1(:) = (eigen1(:)-eigen2(:))/(2.d0*Delta_E)
       eigen2(:) = zI*eigen(:)*conjg(eigen1(:))
       do i=1,max_open
          write(47,9001)abs(eigen(i)),abs(eigen2(i)),i,E
          write(49,9002)real(eigen2(i)),aimag(eigen2(i)),i,E
       end do
       absolute(:)=abs(eigen2(:))
write(40,*)'maxloc',U_copy(:,maxloc(absolute(:),1),nE)
do i=1,max_open
write(34,*)U_copy(i,maxloc(absolute(:),1),nE-1),-U_copy(i,maxloc(absolute(:),1),nE+1)
end do
       write(44,*)absolute(:),nE
       write(43,*)sum(absolute,1),E,maxval(absolute)
       write(42,*)qd(:,nE)
       loc(nE)=maxloc(absolute(:),1) ! ARE YOU SURE IT WORKS?
       write(43,*)'vector=',loc(nE),U(:,loc(nE),nE)

       deallocate(Identity)

       write(45,*)loc(nE),E
       ! Get the eigenphase sum, get resonance parameters

       eigen(:)=log(eigen(:))/(2.d0*zI)/Pi
       eigenphase_sum=0.d0
       eigenphase_sum=sum(real(eigen(:)))
       write(46,*)E,eigenphase_sum
       !--------------------------------------------------
       ! Get dominant partial wave contribution
       max_l=maxloc(U(:,loc(nE),nE)*U(:,loc(nE),nE),1)
       write(41,*)max_l,nE,U(:,loc(nE),nE)*U(:,loc(nE),nE)
       ! sum the l components
       do l=0,lmax
        result =0.d0
        do m=1, 2*l+1
        channel=(l)**2+m
           result= result+(U(channel,loc(nE),nE)*U(channel,loc(nE),nE))
        end do
       write(50,*)result, l, nE
       end do
       !--------------------------------------------------
       deallocate(eigen,eigen1,eigen2,absolute)
       deallocate(S_current)
       deallocate(S1)
       deallocate(S2)
       deallocate(eigen_copy1,eigen_copy2,top1,top2,dot1,dot2)
    end do
9001 format(2f13.5,I,1f13.5)
9002 format(2f13.5,2I)

  end subroutine time_delay2



  subroutine input_matrices(nClosed,max_open,max_index,p_Patch,RHS_co,Z,b,x_coord,y_coord,z_coord)
    use nrtype, only : dpc,dbl,i4b
    use Solve, only : invert_real
    use control, only : option_wf
    implicit none
    integer(kind=i4b),intent(IN) :: nClosed,max_open,max_index,p_Patch
    integer(kind=i4b) :: f,q,p
    real(kind=dbl) :: RHS_co(nClosed,max_open),Z(max_open,max_open),b(max_open),x_coord(nClosed),y_coord(nClosed),z_coord(nClosed)
    real(kind=dbl),dimension(:,:),allocatable :: Z_sav,cc
    write(6,*)'before R matrix input',max_open,nClosed, max_index,(nClosed-max_index)/8,p_Patch

    ! Input Z
RHS_co(:,:)=0.d0
x_coord(:)=0.d0
y_coord(:)=0.d0
z_coord(:)=0.d0

    do f=1,max_open
       read(64,*)b(f)
       !write(74,*)b(f)
       do q=1,max_open
          read(63,*)Z(f,q)
          !write(73,*)Z(f,q)
       end do
    end do
    write(6,*)'after R matrix input',max_open,nClosed, max_index,(nClosed-max_index)/8
    ! Input b

    ! Input RHS_co and coords
if (option_wf=='wf') then
    f=1
    q=0
    do f=1,max_open
       q=0
       do p=1,nClosed-max_index,8
          q=q+1
          read(9000+f,902)x_coord(q),y_coord(q),z_coord(q),RHS_co(p,f)
       end do
    end do
end if
902 format(4e13.5)
    write(6,*)'after cc input',max_open,nClosed
    allocate(Z_sav(max_open,max_open))
    Z_sav=Z
    call invert_real(Z_sav,max_open,max_open)
    RHS_co=RHS_co/sqrt(2.d0)
    !RHS_co=-matmul(RHS_co,Z_sav)
    ! Call Blas 3 subroutine to perform this multiplication
    !RHS_co(:,:)=0.d0
allocate(cc(nClosed,max_open))
cc(:,:)=0.d0
if (option_wf.eq.'wf') call dgemm ( 'N', 'N', p_Patch-max_index, max_open, max_open, -1.d0, RHS_co,p_Patch-max_index, Z, max_open, &
         &                   0.d0, cc, p_Patch-max_index )
RHS_co(:,:)=cc(:,:)
deallocate(cc)
    deallocate(Z_sav)
write(6,*)'end input matrices'
  end subroutine input_matrices


  subroutine  differential(nEmax,max_open,lds,Emax,Emin,loc,max_open_array,l_max,U,qd,S)
    use nrtype, only : dpc,dbl,i4b,Pi
    use Solve
    use gensub
    implicit none
    integer(kind=i4b) :: nEmax,max_open,lds,nE,lmax,max_open_top,i,k
    integer(kind=i4b),dimension(nEmax) :: loc
    real(kind=dbl) :: Emax,Emin,Delta_E,E,eigenphase_sum,kappa,angle,pleg,&
         & differential_cross_section,plgndr
    integer(kind=i4b),dimension(nEmax) :: max_open_array,l_max
    integer(kind=i4b) :: m,mm,mmm,m_0,l,ll,lll,l_0,capL_0,capL,capLL,capLLL,alf
    real(kind=dbl),dimension(max_open,max_open,nEmax) :: U
    real(kind=dbl),dimension(max_open,nEmax) :: qd
    real(kind=dbl),dimension(:),allocatable :: dot1,dot2,AA
    integer(kind=i4b),dimension(:),allocatable :: top1,top2
    complex(kind=dpc),dimension(max_open,max_open,nEmax) :: S
    complex(kind=dpc),dimension(:,:),allocatable :: S_current,S1,S2,Identity
    complex(kind=dpc),dimension(:),allocatable :: eigen,eigen1,eigen2
    complex(kind=dpc),dimension(:),allocatable :: eigen_copy1,eigen_copy2
    complex(kind=dpc), parameter :: zI=(0.d0,1.d0)
    complex(kind=dpc):: effe,coefficient,coeff=0.d0,coeff2=0.d0
    integer(kind=i4b) :: j,index_chan
    integer(kind=i4b), parameter :: max_angle=30
    real(kind=dbl),dimension(:),allocatable :: absolute
    external :: pleg,plgndr
    ! Get differential cross sections
    Delta_E=(Emax-Emin)/dble(nEmax)
    write(6,*)'inside differential',nEmax,max_open
    write(6,*)max_open_array(:),l_max(:)
    max_open_top=max_open_array(nEmax)
    loc(1)=1
    loc(nEmax)=1
    do nE=1,nEmax
       max_open=max_open_array(nE)
       lmax=l_max(nE)
       E=Emin+nE*Delta_E
       kappa=sqrt(2.d0*E)
       !write(6,*)max_open,lmax,E,Emin,Emax,Delta_E
       allocate(Identity(max_open_top,max_open_top))
       allocate(eigen(max_open_top),eigen1(max_open_top),eigen2(max_open_top),absolute(max_open_top))
       allocate(S_current(max_open_top,max_open_top))
       allocate(S1(max_open_top,max_open_top))
       allocate(S2(max_open_top,max_open_top))

       Identity(:,:)=(0.d0,0.d0)
       do j=1,max_open_top
          Identity(j,j)=(1.d0,0.d0)
       end do
       S_current(1:max_open,1:max_open)=S(1:max_open,1:max_open,nE)
       !S1(:,:)=S(:,:,nE-1)
       !S2(:,:)=S(:,:,nE+1)
       !eigen(:)=exp(2.d0*zI*Pi*qd(:,nE))
       do i=1,max_open
         eigen(i)=exp(2.d0*zI*Pi*qd(i,nE))
       end do

       !do i=max_open+1, max_open_top
       !   S_current(i,i)=cmplx(1.d0,0.d0)
       !   eigen(i)=cmplx(1.d0,0.d0)
       !end do

       !-----------------------------------------
       ! Simple
       !do i=1,max_angle
       !angle=(dble(i-1)*(Pi/dble(max_angle)))
       !effe=0.d0
       !do j=0,l_max(nEmax)-1
       !effe=effe+(eigen(j**2+2*j+1)-1)*(2*j+1)/(2.d0*zI*kappa)*plgndr(j,0,cos(angle))
       !end do
       !differential_cross_section=effe*conjg(effe)
       !write(99+nE,*)angle*180.d0/Pi,differential_cross_section
       !end do
       ! Complicated
       allocate(AA(0:lmax))
       Identity(:,:)=cmplx(0.d0,0.d0)
       do j=1,max_open_top
          Identity(j,j)=eigen(j)
       end do
       !S_current=matmul(transpose(U(:,:,nE)),Identity)
       !S_current=matmul(S_current,U(:,:,nE))
       do i=1,max_open
          S_current(i,i)=1.d0-S_current(i,i)
          !write(29,*)S_current(i,i),nE
       end do

       alf=2
       !---------------------------------------------
       ! Initialization of Clebsch-Gordan subroutines
       call setup()
       !---------------------------------------------
       do k=0,lmax
          write(6,*)'l=',k
          coeff=0.d0
          coefficient=0.d0
          do capL_0=1,max_open
             l_0=int(sqrt(capL_0-1.d0)); m_0=capL_0-l_0*(l_0+1)-1
             write(28,*)capL_0,l_0,m_0
             do capL=1,max_open
                l=int(sqrt(capL-1.d0)); m=capL-l*(l+1)-1
                coeff2=0.d0
                do capLL=1,max_open
                   ll=int(sqrt(capLL-1.d0)); mm=capLL-ll*(ll+1)-1
                   coeff=0.d0
                   do capLLL=1,max_open
                      lll=int(sqrt(capLLL-1.d0)); mmm=capLLL-lll*(lll+1)-1
                      !          coefficient=coefficient+zI**(l_0-l)*(-1.d0)**m_0*sqrt((2*l_0+1.d0)*(2*l+1.d0))*(S_current(capL_0,capL))*&
                      !            & zI**(lll-ll)*(-1.d0)**m*sqrt((2*ll+1.d0)*(2*lll+1.d0))*(conjg(S_current(capLL,capLLL)))*&
                      !            & clebsh(alf*l,alf*lll,alf*k,0,0,0)*clebsh(alf*l_0,alf*ll,alf*k,0,0,0)*&
                      !            & clebsh(alf*l_0,alf*ll,alf*k,-alf*m_0,alf*mm,alf*(mm-m_0))*&
                      !            &clebsh(alf*l,alf*lll,alf*k,-alf*m,alf*mmm,alf*(mm-m_0))
                      coeff=coeff+ &
                           & zI**(lll-ll)*sqrt((2*lll+1.d0))*(conjg(S_current(capLL,capLLL)))*&
                           & clebsh(alf*l,alf*lll,alf*k,0,0,0)*&
                           & clebsh(alf*l,alf*lll,alf*k,-alf*m,alf*mmm,alf*(mm-m_0))*(-1.d0)**(mmm)
                   end do
                   coeff2=coeff2+coeff*clebsh(alf*l_0,alf*ll,alf*k,-alf*m_0,alf*mm,alf*(mm-m_0))*&
                        & clebsh(alf*l_0,alf*ll,alf*k,0,0,0)*sqrt((2*ll+1.d0))
                end do
                coefficient=coefficient+coeff2*&
                     & zI**(l_0-l)*(-1.d0)**(m_0)*&
                     & sqrt((2*l_0+1.d0)*(2*l+1.d0))*(S_current(capL_0,capL))
             end do
          end do
          AA(k)=real(coefficient)/(2.d0*k+1)
          write(20,*)real(coefficient),aimag(coefficient),k,E
       end do
       write(21,*)Pi/(2.d0*E)*(AA(0)-AA(1)/3.d0),E,Pi/(2.d0*E)*(AA(0))
       do i=1,max_angle+1
          angle=(dble(i-1)*(Pi/dble(max_angle)))
          differential_cross_section=0.d0
          do k=0,lmax
             differential_cross_section=differential_cross_section+AA(k)*plgndr(k,0,cos(angle))/(8.d0*E)
          end do
          write(99+nE,*)angle*180.d0/Pi,differential_cross_section
       end do

       deallocate(AA)
       !-------------------------------------------
       write(6,*)'after allocs'
       deallocate(eigen,eigen1,eigen2,absolute)
       deallocate(S_current)
       deallocate(S1)
       deallocate(S2)
       deallocate(Identity)
    end do
    stop
9001 format(2f13.5,I,1f13.5)
  end subroutine differential


end module alpha
