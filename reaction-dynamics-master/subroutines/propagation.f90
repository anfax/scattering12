!!!*************************************************************
! 文件/File: propagation.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: propagation.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:45
!*************************************************************

! Propogation of the wave
subroutine propagation(psi0,Rgrid,r1grid,thetagrid,tkk1)
use constants
    implicit none
    real*8,dimension(1:nR)::Rgrid,DampR,kRgrid
    real*8,dimension(1:nr1)::r1grid,Dampr1,kr1grid
    real*8,dimension(1:ntheta)::thetagrid
    real*8,dimension(1:ntheta,1:ntheta)::tkk1
    real*8,dimension(1:nR,1:nr1,1:ntheta):: Tgrid,psipot,psi0,psiang,psirad,psitot,psif,psik
    real*8,dimension(1:nR,1:nr1):: damp

    real*8,dimension(3)::Rad !The distances between the atoms,consists of r1,r2,r3 in anticlockwise order(atoms a,b,c)
    real*8,dimension(1:nR,1:nr1,1:ntheta)::vpot !The Potential
    real*8 Jtot,js,Kcap,total,vpott

    integer*8 i,j,k,nc,k1
    character(20)::filename


    Jtot = 1.0
    Kcap = 0.0
    !print*,nR,Rgrid

    ! Calculating the Tgrid (Angular Kinetic energy Operator)
    do i=1,nR
        do j=1,nr1
            do k=1,ntheta
                js = k-1
                    total = (hcross**2/2.0) *(( (Jtot*(Jtot+1) + js*(js+1) - 2.0*Kcap**2)/(rmR*Rgrid(i)**2)) &
                    + (js*(js+1))/(rmr1*r1grid(j)**2) )
                    if(total .ge. 5.00/27.2116) total = 5.00/27.2116 ! capping the value in Tgrid
                Tgrid(i,j,k) = total
            end do
        end do
    end do
    print*, "(Propogation)-> Tgrid (Angular Kinetic energy Operator) calculated"

    do i=1,nR
        do j=1,nr1
            do k=1,ntheta
                !Initialisation of the grid from Jackobi coordinates  (cmass/(bmass+cmass))
                !Assumption- All three masses are equal
                Rad(1) = sqrt(Rgrid(i)**2 + (r1grid(j)*0.5)**2 + &
                  2.0*r1grid(j)*Rgrid(i)*cos(thetagrid(k))*0.5)
                Rad(2) = r1grid(j)
                Rad(3) = sqrt(Rgrid(i)**2 + (r1grid(j)*0.5)**2 - &
                  2.0*r1grid(j)*Rgrid(i)*cos(thetagrid(k))*0.5)
                !Calling bkmp2 for potential grid
                call bkmp2(Rad,vpott)
                vpot(i,j,k) = vpott + 0.17449577081d0
                if(vpot(i,j,k) .ge. 0.2) vpot(i,j,k) = 0.2
            enddo
        enddo 
    enddo 

    open(unit=11,file='potential.out')  
    do i=1,nR
        do j=1,nr1
            write(11,*)Rgrid(i),r1grid(j),vpot(i,j,20)
        enddo
        write(11,*) 
    enddo 
    close(11)

    print*, "(Propogation)-> Grid Initialisation from Jackobi coordinates done"

    ! The damping function 
    call damping(DampR,Rgrid,nR)
    call damping(Dampr1,r1grid,nr1)

    open(unit=20,file='damping.out')
    do i=1,nR
        do j=1,nr1
            damp(i,j)=DampR(i)*Dampr1(j) 
!           write(20,*)Rgrid(i),r1grid(j),damp(i,j)
        enddo 
!       write(20,*)
    enddo 
    close(20)

    ! Momentum grid points for Rgrid and r1grid
    ! used in calculating Radial kinetic energy
    call makeMomentumGrid(kRgrid,nR,Rmax-Rmin)
    call makeMomentumGrid(kr1grid,nr1,r1max-r1min)
    print*, "(Propogation)-> Momentum grid points for Rgrid and r1grid made"


    !!!!!Loop over chebyshevs (ncheb in modules)
    Do nc=1,ncheb

      ! Potential application
      do i=1,nR
        do j=1,nr1
            do k=1,ntheta
                psipot(i,j,k)= psi0(i,j,k)*vpot(i,j,k) 
            enddo
        enddo 
      enddo 

      ! Going from dvr to fbr and back
      ! Rotational Kinetic energy part of the propogation
      do i=1,nR
        do j=1,nr1
            do k=1,ntheta
                total=0.0
                do k1=1,ntheta
                    total=total+tkk1(k,k1)*psi0(i,j,k)
                enddo
                psiang(i,j,k)=total
            enddo
        enddo 
      enddo 

! Operator
      do i=1,nR
        do j=1,nr1
            do k=1,ntheta
                  psiang(i,j,k) = Tgrid(i,j,k) * psiang(i,j,k)
            end do
        end do
      end do

      do i=1,nR
        do j=1,nr1
            do k=1,ntheta
                total=0.0
                do k1=1,ntheta
                    total=total+tkk1(k1,k)*psiang(i,j,k) ! Transpose of the original tkk1 matrix
                enddo
                psiang(i,j,k)=total
            enddo
        enddo 
      enddo 

      call kineticRad(psi0,krgrid,kr1grid,psirad)

      
      ! H -> Hscaled
      !psitot = (psiang + psipot + psirad - (Hmax+Hmin)/2.0 )/((Hmax-Hmin)/2.0)
      psitot = ( psipot - (Hmax+Hmin)/2.0 )/((Hmax-Hmin)/2.0)

! Potential problem (need to diagonalise for finding psi)

      ! Recursively calling chebyshevs after this
      if (nc .eq. 1) then
        do i=1,nR
            do j=1,nr1
                do k=1,ntheta
                    psif(i,j,k)=psitot(i,j,k)!*damp(i,j) 
                enddo
            enddo  
        enddo 

      else 

        do i=1,nR
            do j=1,nr1
                do k=1,ntheta
                    psif(i,j,k)= ( 2*psitot(i,j,k) - psik(i,j,k)*damp(i,j) )*damp(i,j) 
                enddo
            enddo 
        enddo 

      endif

      psik = psi0
      psi0 = psif

!      if (mod(nc,100) .eq. 0) then
          write(filename,7)nc
          7 format('data/wav',i0,'.out')
          open(unit=3,file=filename)
          do i=1,nR
              do j=1,nr1
                write(3,*),Rgrid(i),r1grid(j),psi0(i,j,1)**2
              enddo
          write(3,*)
          end do
          close(3)
!      endif

    Enddo

end subroutine
