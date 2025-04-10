program triatomic
use constants
use angulardvr
use legendre

    implicit none
    integer *8 i,j,k
    real*8,dimension(:),allocatable::Rgrid,r1grid,thetagrid,weight,evalr1,Rfunc,bcwav
    real*8,dimension(:,:),allocatable::evecr1,tkk1
    real*8,dimension(:,:,:),allocatable::psi
    real*8::NormR,Normr1

    allocate(Rgrid(1:nR))
    allocate(Rfunc(1:nR))
    allocate(r1grid(1:nr1))
    allocate(bcwav(1:nr1))
    allocate(thetagrid(1:ntheta))
    allocate(weight(1:ntheta))
    allocate(evalr1(1:nr1))
    allocate(evecr1(1:nr1,1:nr1))
    allocate(psi(1:nr,1:nr1,1:ntheta))
    allocate(tkk1(1:ntheta,1:ntheta))

    open(unit=9,file='gauss.out')
    do i=1,nR
        Rgrid(i)=Rmin+(i-1)*delR
    end do

    NormR = 0.0
    do i=1,nR
        Rfunc(i)=exp(-((Rgrid(i)-Ro)**2)/(2*(dell**2)))*cos(ko*Rgrid(i))
        NormR = NormR + Rfunc(i)**2
        write(9,*),Rgrid(i),Rfunc(i)**2
    end do
    print*, "Normalisation R :",NormR

    do i=1,nr1
        r1grid(i)=r1min+(i-1)*delr1
    end do

    call bcpsi(nr1,r1grid,delr1,evalr1,evecr1) !psi of vibration
    
    Normr1 = 0.0
    open (unit=12,file='bcpsi.out')
    do i=1,nr1
        bcwav(i) = evecr1(i,vstate+1)
        Normr1 = Normr1 + bcwav(i)**2
        write(12,*)r1grid(i),bcwav(i)
    end do
    print*,"Normalisation r1: ",Normr1

    call angle(thetagrid,weight,mstate) !grid of angle

    open (unit=5,file='thetagrid.out')
    do j=1,ntheta
        write(5,*)thetagrid(j),Plgndr(lstate,mstate,cos(thetagrid(j))),weight(j)
    end do

    open(unit=23,file='wav0.out')
     do i=1,nR
      do j=1,nr1
            do k=1,ntheta
              psi(i,j,k)=Rfunc(i)*bcwav(j)*weight(k)*Plgndr(lstate,mstate,cos(thetagrid(k)))/sqrt(Normr1*NormR)
            enddo
        end do
     enddo

      do i=1,nR
      do j=1,nr1
         write(23,*),Rgrid(i),r1grid(j),psi(i,j,1)**2
        enddo
         write(23,*)
     end do

     ! Making the tkk1 matrix from weights, to be used in propogation
     do j=0,ntheta-1 !j= lstate
        do i=1,ntheta ! i is β 
            tkk1(j+1,i) = sqrt(weight(i)) * Plgndr(j+1,mstate,cos(thetagrid(i))) ! assuming that Ω is mstate
        end do
    end do

    call propagation(psi,Rgrid,r1grid,thetagrid,tkk1)

end program triatomic
