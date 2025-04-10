!!!*************************************************************
! 文件/File: dbs3ine90.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: dbs3ine90.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

subroutine dbs3ines(xdata,ydata,zdata,fdata,nxdata,nydata,nzdata,nxknot,nyknot,nzknot,xvec,yvec,zvec,nvec,valore)

use numeric

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! IMSL name:  dbs3in (double precision version)
!
! purpose:    compute a three-dimensional tensor-product spline
!             interpolant, returning the tensor-product B-spline
!             coefficients.
!
! usage:      call dbs3in(nxdata, xdata, nydata, ydata, nzdata,
!                         zdata, fdata, ldf, mdf, kxord, kyord,
!                         kzord, xknot, yknot, zknot, bscoef)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! List of variables in the call
! x-y-zdata : grid points
! nx-y-zdata : dimension of the grid

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      use bspline
 
      implicit none

!
!        specifications for parameters
!

      integer    kxord, kyord, kzord, ldf, mdf, nxdata, nxknot, nxvec,  &
     &           nydata, nyknot, nyvec, nzdata, nzknot, nzvec,order,nvec

      integer    i, j, k, nxcoef, nycoef, nzcoef
!      double precision  bscoef(nxdata,nydata,nzdata), f,                &
!     &           fdata(nxdata,nydata,nzdata), &
!     &           , x, xdata(nxdata), xknot(nxknot), xvec(nvec), y,     &
!     &           ydata(nydata), yknot(nyknot), yvec(nvec), z,          &
!     &           zdata(nzdata), zknot(nzknot), zvec(nvec),valore(nvec)
      real(kind=dabl) :: val,f,x,y,z
      real(kind=dabl), dimension(nxdata),intent(in) :: xdata
      real(kind=dabl), dimension(nydata),intent(in) :: ydata 
      real(kind=dabl), dimension(nzdata),intent(in) :: zdata
      real(kind=dabl), dimension(nvec),intent(in) :: xvec
      real(kind=dabl), dimension(nvec),intent(in) :: yvec
      real(kind=dabl), dimension(nvec),intent(in) :: zvec
      real(kind=dabl), dimension(nxdata,nydata,nzdata),intent(in) :: fdata
      real(kind=dabl), dimension(nvec),intent(out) :: valore
      real(kind=dabl), dimension(:,:,:),allocatable :: bscoef
      real(kind=dabl), dimension(:),allocatable :: xknot
      real(kind=dabl), dimension(:),allocatable :: yknot
      real(kind=dabl), dimension(:),allocatable :: zknot

!      real(kind=dabl) :: dbs3vl
!      external :: dbs3vl  
!
!        define function.
!
!        define dimensions of fdata

      ldf =nxdata
      mdf=nydata
      order=nxknot-nxdata
      write(6,*)order
      write(33,*)xvec
      kxord=order
      kyord=order
      kzord=order
      nxvec=nvec
      nyvec=nvec
      nzvec=nvec
!
!        generate knots
!

allocate(xknot(nxknot))
allocate(yknot(nyknot))
allocate(zknot(nzknot))

write(36,*)xdata
      call dbsnak (nxdata, xdata, kxord, xknot)
      call dbsnak (nydata, ydata, kyord, yknot)
      call dbsnak (nzdata, zdata, kzord, zknot)

!
!        interpolate
!
!do i=1,nxdata
!do j=1,nydata
!do k=1,nzdata
!alpha=fdata(i,j,k)
!beta=bscoef(i,j,k)
!end do 
!end do
!end do
write(6,*)'before allocate bscoef'
allocate(bscoef(nxdata,nydata,nzdata))
write(6,*)'after allocate bscoef'

      call dbs3in (nxdata, xdata, nydata, ydata, nzdata, zdata, fdata,  &
     &            ldf, mdf, kxord, kyord, kzord, xknot, yknot, zknot,   &
     &            bscoef)
write(6,*)'after dbs3in'

      nxcoef = nxdata
      nycoef = nydata
      nzcoef = nzdata

!
!        write heading
!

      write (6,99999)

!
!        call the evaluation routine.
!

!      call dbs3gd (0, 0, 0, nxvec, xvec, nyvec, yvec, nzvec, zvec,      &
!     &            kxord, kyord, kzord, xknot, yknot, zknot, nxcoef,     &
!     &            nycoef, nzcoef, bscoef, value, nxvec, nyvec)
do i=1,nvec
 x=xvec(i)
 y=yvec(i)
 z=zvec(i)
  val=dbs3vl(x,y,z,kxord,kyord,kzord,xknot,yknot,zknot,nxdata,nydata,nzdata,bscoef)
valore(i)=val
end do

deallocate(bscoef)
deallocate(xknot)
deallocate(yknot)
deallocate(zknot)

      do i = 1, nxvec
         do j = 1, nyvec
            do k = 1, nzvec
!               write (6,'(4f13.4, f13.6)') xvec(i), yvec(k) ,            &
!     &                                       zvec(k), value(i,j,k),     &
!     &                                       f(xvec(i),yvec(j),zvec(k)) &
!     &                                        , value(i,j,k)
            end do
         end do
      end do

99999 format (10x, 'x', 11x, 'y', 10x, 'z', 10x, 's(x,y,z)', 7x,        &
     &       'error')

!     Evaluate the function at the point requested using De Boor's subroutines

!ndimx=nxdata+kxord
!ndimy=nydata+kyord
!ndimz=nzdata+kzord

!do element=1,n

!  call Spline_value(ndimx,ndimy,ndimz,nxdata,nydata,nzdata,xdata,ydata,zdata,bscoef,xtx,xty,xtz,order)
!end do
      end subroutine dbs3ines

