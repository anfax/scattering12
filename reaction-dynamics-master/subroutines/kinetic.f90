!!!*************************************************************
! 文件/File: kinetic.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: kinetic.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:45
!*************************************************************

subroutine kineticRad(psi0,kRgrid,kr1grid,psiout)
    use constants
    implicit none
    real*8,dimension(1:nR,1:nr1,1:ntheta)::psi0,psiR,psir1,psiout
    real*8,dimension(1:nR)::kRgrid
    real*8,dimension(1:nr1)::kr1grid
    integer*8 info,i,j,k
    real*8,dimension(:),allocatable::comm,temp

    info = 0

    do i=1,nR
     do j=1,nr1
      do k=1,ntheta
        psiR(i,j,k) = psi0(i,j,k)
        psir1(i,j,k) = psi0(i,j,k)
      enddo
     enddo 
    enddo 

    ! --------------- For R grid ----------------
    allocate(comm(1:3*nR+300))
    allocate(temp(nR))

    ! Changing to momentum coordinates
    do k=1,ntheta
        do j=1,nr1
            do i=1,nR
                temp(i) = psiR(i,j,k)
            end do
            call DZFFT(2,nR,temp,comm,info)
            if(info.lt.0) print*,"(DZFFT-R)-> Error code ",info
            do i=1,nR
                psiR(i,j,k) = temp(i)
            end do
        end do
    end do

    ! Kinetic Energy calculation
    do i=1,nR
        do j=1,nr1
            do k=1,ntheta
                psiR(i,j,k) = psiR(i,j,k) * (-1*kRgrid(i)**2/(2.0*rmR) )
            end do
        end do
    end do

    ! Back from momentum coordinates
    do k=1,ntheta
        do j=1,nr1
            do i=1,nR
                temp(i) = psiR(i,j,k)
            end do
            call ZDFFT(2,nR,temp,comm,info)
            if(info.lt.0) print*,"(ZDFFT-R)-> Error code ",info
            do i=1,nR
                psiR(i,j,k) = temp(i)
            end do
        end do
    end do
    deallocate(comm)
    deallocate(temp)


    ! --------------- For r1 grid ----------------
    allocate(comm(1:3*nr1+300))
    allocate(temp(nr1))
    
    ! Changing to momentum coordinates
    do k=1,ntheta
        do i=1,nR
            do j=1,nr1
                temp(j) = psir1(i,j,k)
            end do
            call DZFFT(2,nr1,temp,comm,info)
            if(info.lt.0) print*,"(DZFFT-r)-> Error code ",info
            do j=1,nr1
                psir1(i,j,k) = temp(j)
            end do
        end do
    end do

    ! Kinetic Energy calculation
    do i=1,nR
        do j=1,nr1
            do k=1,ntheta
                psir1(i,j,k) = psir1(i,j,k) * ( -1*kr1grid(i)**2/(2.0*rmr1) )
            end do
        end do
    end do

    ! Back from momentum coordinates
    do k=1,ntheta
        do i=1,nR
            do j=1,nr1
                temp(j) = psir1(i,j,k)
            end do
            call ZDFFT(2,nr1,temp,comm,info)
            if(info.lt.0) print*,"(ZDFFT-r)-> Error code ",info
            do j=1,nr1
                psir1(i,j,k) = temp(j)
            end do
        end do
    end do
    deallocate(comm)
    deallocate(temp)

    do i=1,nR
     do j=1,nr1
      do k=1,ntheta
       psiout(i,j,k) = psiR(i,j,k) + psir1(i,j,k) 
      enddo
     enddo 
    enddo 

end subroutine kineticRad