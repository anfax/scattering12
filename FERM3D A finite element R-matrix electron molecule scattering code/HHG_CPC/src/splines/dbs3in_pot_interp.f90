!!!*************************************************************
! 文件/File: dbs3in_pot_interp.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: dbs3in_pot_interp.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

module potential_interp
use numeric,only : dbl=>dabl,i4b=>i4
      integer(kind=i4b) :: nnxdata,nnydata,nnzdata
end module potential_interp

subroutine dbs3ines(xdata,ydata,zdata,fdata,nxdata,nydata,nzdata,nxknot,nyknot,nzknot,bscoef,xknot,yknot,zknot)
use potential_interp
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

      integer(kind=i4)    kxord, kyord, kzord, ldf, mdf, nxdata, nxknot, nxvec,  &
     &           nydata, nyknot, nyvec, nzdata, nzknot, nzvec,order,nvec

      integer(kind=i4)    i, j, k, nxcoef, nycoef, nzcoef
!      double precision  bscoef(nxdata,nydata,nzdata), f,                &
!     &           fdata(nxdata,nydata,nzdata), &
!     &           , x, xdata(nxdata), xknot(nxknot), xvec(nvec), y,     &
!     &           ydata(nydata), yknot(nyknot), yvec(nvec), z,          &
!     &           zdata(nzdata), zknot(nzknot), zvec(nvec),valore(nvec)
      real(kind=dabl) :: val,f,x,y,z
      real(kind=dabl), dimension(nxdata),intent(in) :: xdata
      real(kind=dabl), dimension(nydata),intent(in) :: ydata 
      real(kind=dabl), dimension(nzdata),intent(in) :: zdata
!      real(kind=dabl), dimension(nvec),intent(in) :: xvec
!      real(kind=dabl), dimension(nvec),intent(in) :: yvec
!      real(kind=dabl), dimension(nvec),intent(in) :: zvec
      real(kind=dabl), dimension(nxdata,nydata,nzdata), intent(in) :: fdata
      real(kind=dabl), dimension(nxdata,nydata,nzdata), intent(out) :: bscoef
      real(kind=dabl), dimension(nxknot),intent(out) :: xknot
      real(kind=dabl), dimension(nyknot),intent(out) :: yknot
      real(kind=dabl), dimension(nzknot),intent(out) :: zknot
!      real(kind=dabl) :: dbs3vl
!      external :: dbs3vl  
!
!        define function.
!
!        define dimensions of fdata
nnxdata=nxdata
nnydata=nydata
nnzdata=nzdata

      ldf =nxdata
      mdf=nydata
      order=nxknot-nxdata
      !write(6,*)order
      kxord=order
      kyord=order
      kzord=order
      nxvec=nvec
      nyvec=nvec
      nzvec=nvec
!
!        generate knots
!

!write(36,*)nxdata,xdata,nydata,ydata,nzdata,zdata
      call dbsnak (nxdata, xdata, kxord, xknot)
!write(39,*)nxknot,xknot
      call dbsnak (nydata, ydata, kyord, yknot)
!write(40,*)nyknot,yknot
      call dbsnak (nzdata, zdata, kzord, zknot)
!write(41,*)nzknot,zknot
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
!write(6,*)'before dbs3in',nxdata, xdata, nydata, ydata, nzdata, zdata,ldf, mdf, kxord, kyord, kzord

      call dbs3in (nxdata, xdata, nydata, ydata, nzdata, zdata, fdata,  &
     &            ldf, mdf, kxord, kyord, kzord, xknot, yknot, zknot,   &
     &            bscoef)
!write(6,*)'after dbs3in'

      nxcoef = nxdata
      nycoef = nydata
      nzcoef = nzdata

!
!        write heading
!

!      write (6,99999)

!
!        call the evaluation routine.
!

!      call dbs3gd (0, 0, 0, nxvec, xvec, nyvec, yvec, nzvec, zvec,      &
!     &            kxord, kyord, kzord, xknot, yknot, zknot, nxcoef,     &
!     &            nycoef, nzcoef, bscoef, value, nxvec, nyvec)
!do i=1,nvec
! x=xvec(i)
! y=yvec(i)
! z=zvec(i)
!  val=dbs3vl(x,y,z,kxord,kyord,kzord,xknot,yknot,zknot,nxdata,nydata,nzdata,bscoef)
!valore(i)=val
!end do


!XXXXXXXXXXXXXXXXXXX
!deallocate(bscoef)
!deallocate(xknot)
!deallocate(yknot)
!deallocate(zknot)
!XXXXXXXXXXXXXXXXXXXX

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

