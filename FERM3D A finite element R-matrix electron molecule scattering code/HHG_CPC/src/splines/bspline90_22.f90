!!!*************************************************************
! 文件/File: bspline90_22.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: bspline90_22.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************

module numeric 

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!   VERSION 2.2
!
!   f90 VERSION
!
!   This library contains routines for B-spline interpolation in
!   one, two, and three dimensions. Part of the routines are based
!   on the book by Carl de Boor: A practical guide to Splines (Springer,
!   New-York 1978) and have the same calling sequence and names as
!   the corresponding routines from the IMSL library. For documen-
!   tation see the additional files. NOTE: The results in the demo
!   routines may vary slightly on different architectures.
!
!   by W. Schadow 12/04/99
!   last changed by W. Schadow 07/28/2000
!
!
!   Wolfgang Schadow
!   TRIUMF
!   4004 Wesbrook Mall
!   Vancouver, B.C. V6T 2A3
!   Canada
!
!   email: schadow@triumf.ca  or  schadow@physik.uni-bonn.de
!
!   www  : http://www.triumf.ca/people/schadow
!
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!   Copyright (C) 2000 Wolfgang Schadow
!
!   This library is free software; you can redistribute it and/or
!   modify it under the terms of the GNU Library General Public
!   License as published by the Free Software Foundation; either
!   version 2 of the License, or (at your option) any later version.
!
!   This library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!   Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public
!   License along with this library; if not, write to the
!   Free Software Foundation, Inc., 59 Temple Place - Suite 330,
!   Boston, MA  02111-1307, USA.
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  integer, parameter :: sgl = kind(1.0)
  integer, parameter :: dabl = kind(1.0d0)
  INTEGER, PARAMETER :: I4 = SELECTED_INT_KIND(9)
end module numeric


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


module bspline

!
!  ------------------------------------------------------------------
!
!
!   The following routines are included:
!
!            dbsnak
!
!            dbsint
!            dbsval
!            dbsder
!            dbs1gd
!
!            dbs2in
!            dbs2dr
!            dbs2vl
!            dbs2gd
!
!            dbs3in
!            dbs3vl
!            dbs3dr
!            dbs3gd
!
!  ------------------------------------------------------------------
!

  private

  public dbsnak
  public dbsint, dbsval, dbsder, dbs1gd
  public dbs2in, dbs2dr, dbs2vl, dbs2gd
  public dbs3in, dbs3vl, dbs3dr, dbs3gd


contains


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbsnak(nx,xvec,kxord,xknot)

!
!  Compute the `not-a-knot' spline knot sequence.
!  (see de Boor p. 167)
!
!   nx     - number of data points.  (input)
!   xvec   - array of length ndata containing the location of the
!            data points.  (input)
!   kxord  - order of the spline.  (input)
!   xknot  - array of length ndata+korder containing the knot
!            sequence.  (output)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in) :: nx, kxord

    real(kind=dabl), dimension(nx), intent(in)        :: xvec
    real(kind=dabl), dimension(nx+kxord), intent(out) :: xknot

    real(kind=dabl) :: eps
    integer(kind=i4)        :: ix
    logical        :: first = .true.

    save first,eps

    if (first) then
       first=.false.
       eps = epsilon(1.0_dabl)
       write(6,*) "subroutine dbsnak: "
       write(6,*) "eps = ",eps
    endif

    if((kxord .lt. 0) .or. (kxord .gt. nx)) then
       write(6,*) "subroutine dbsnak: error"
       write(6,*) "0 <= kxord <= nx is required."
       write(6,*) "kxord = ", kxord, " and nx = ", nx,  " is given."
       stop
    endif

    do ix = 1, kxord
       xknot(ix) = xvec(1)
    end do
    if(mod(kxord,2) .eq. 0) then
       do ix = kxord+1, nx
          xknot(ix) = xvec(ix-kxord/2)
       end do
    else
       do ix = kxord+1, nx
          xknot(ix) = 0.5_dabl * (xvec(ix-kxord/2) + xvec(ix-kxord/2-1))
       end do
    endif
    do ix = nx+1, nx+kxord
       xknot(ix) = xvec(nx) * (1.0_dabl + eps)
    end do
  end subroutine dbsnak


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbsint(nx,xvec,xdata,kx,xknot,bcoef)

!
!  Computes the spline interpolant, returning the B-spline coefficients.
!  (see de Boor p. 204)
!
!   nx     - number of data points.  (input)
!   xvec   - array of length nx containing the data point
!            abscissas.  (input)
!   xdata  - array of length ndata containing the data point
!            ordinates.  (input)
!   kx     - order of the spline.  (input)
!            korder must be less than or equal to ndata.
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   bscoef - array of length ndata containing the B-spline
!            coefficients.  (output)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                          :: nx, kx
    real(kind=dabl), dimension(nx), intent(in)    :: xdata, xvec
    real(kind=dabl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dabl), dimension(nx), intent(out)   :: bcoef

    integer(kind=i4)                                :: nxp1, kxm1, kpkm2, leftx, lenq
    integer(kind=i4)                                :: ix, ik,ilp1mx, jj, iflag
    real(kind=dabl)                         :: xveci
    real(kind=dabl), dimension(:),allocatable :: work
allocate(work((2*kx-1)*nx))

    nxp1  = nx + 1
    kxm1  = kx - 1
    kpkm2 = 2 * kxm1
    leftx = kx
    lenq  = nx * (kx + kxm1)

    do ix = 1, lenq
       work(ix) = 0.0_dabl
    end do

    do  ix = 1, nx
       xveci  = xvec(ix)
       ilp1mx = min0(ix+kx,nxp1)
       leftx   = max0(leftx,ix)
       if (xveci .lt. xknot(leftx)) goto 998
30     if (xveci .lt. xknot(leftx+1)) go to 40
       leftx = leftx + 1
       if (leftx .lt. ilp1mx) go to 30
       leftx = leftx - 1
       if (xveci .gt. xknot(leftx+1)) goto 998
40     call bsplvb (xknot,nx+kx,kx,1,xveci,leftx,bcoef)
       jj = ix - leftx + 1 + (leftx - kx) * (kx + kxm1)
       do ik = 1, kx
          jj       = jj + kpkm2
          work(jj) = bcoef(ik)
       end do
    end do

    call banfac(work,kx+kxm1,nx,kxm1,kxm1,iflag)

    if (iflag .ne. 1) then
       write(6,*) "subroutine dbsint: error"
       write(6,*) "no solution of linear equation system !!!"
       stop
    end if

    do ix = 1, nx
       bcoef(ix) = xdata(ix)
    end do

    call banslv(work,kx+kxm1,nx,kxm1,kxm1,bcoef)
deallocate(work)
    return

998 write(6,*) "subroutine dbsint:"
    write(6,*) "xknot(ix) <= xknot(ix+1) required."
    write(6,*) ix,xknot(ix),xknot(ix+1)

    stop

  end subroutine dbsint


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbsval(x,kx,xknot,nx,bcoef)

!
!  Evaluates a spline, given its B-spline representation.
!
!   x      - point at which the spline is to be evaluated.  (input)
!   kx     - order of the spline.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   nx     - number of B-spline coefficients.  (input)
!   bcoef  - array of length nx containing the B-spline
!            coefficients.  (input)
!   dbsval - value of the spline at x.  (output)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                          :: nx, kx
    real(kind=dabl)                               :: dbsval
    real(kind=dabl)                               :: x
    real(kind=dabl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dabl), dimension(nx), intent(in)    :: bcoef

    integer(kind=i4)                       :: il, ik, ix, leftx
    real(kind=dabl)                :: save1, save2
    real(kind=dabl), dimension(:),allocatable :: work, dl, dr

allocate(work(kx),dl(kx),dr(kx))
!
!     check if xknot(i) <= xknot(i+1) and calculation of i so that
!     xknot(i) <= x < xknot(i+1)
!

    leftx = 0

    do ix = 1,nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(6,*) "subroutine dbsval:"
          write(6,*) "xknot(ix) <= xknot(ix+1) required."
          write(6,*) ix,xknot(ix),xknot(ix+1)
          stop
       endif
       if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
    end do

    if(leftx .eq. 0) then
       write(6,*) "subroutine dbsval:"
       write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
       write(6,*) "x = ", x
       stop
    endif

    do ik = 1, kx-1
       work(ik) = bcoef(leftx+ik-kx)
       dl(ik)   = x - xknot(leftx+ik-kx)
       dr(ik)   = xknot(leftx+ik) - x
    end do

    work(kx)  = bcoef(leftx)
    dl(kx)    = x - xknot(leftx)

    do ik = 1, kx-1
       save2 = work(ik)
       do il = ik+1, kx
          save1 = work(il)
          work(il) = (dl(il) * work(il) + dr(il-ik) * save2)                  &
               &           / (dl(il) + dr(il - ik))
          save2 = save1
       end do
    end do

    dbsval = work(kx)
deallocate(work,dl,dr)

  end function dbsval


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbsder(iderx,x,kx,xknot,nx,bcoef)

!
!  Evaluates the derivative of a spline, given its B-spline representation.
!
!
!   iderx  - order of the derivative to be evaluated.  (input)
!            in particular, iderx = 0 returns the value of the
!            spline.
!   x      - point at which the spline is to be evaluated.  (input)
!   kx     - order of the spline.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   nx     - number of B-spline coefficients.  (input)
!   bcoef  - array of length nx containing the B-spline
!            coefficients.  (input)
!   dbsder - value of the iderx-th derivative of the spline at x.
!            (output)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                          :: iderx, kx, nx
    real(kind=dabl)                               :: dbsder
    real(kind=dabl), intent(in)                   :: x
    real(kind=dabl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dabl), dimension(nx), intent(in)    :: bcoef

    integer(kind=i4)                       :: ix, ik, il, leftx
    real(kind=dabl)                :: save, save1, save2, y, sum, dik
    real(kind=dabl), dimension(:),allocatable :: work, dl, dr,bsp
allocate(work(kx),dl(kx),dr(kx),bsp(kx))

!
!     check if xknot(i) <= xknot(i+1) and calculation of i so that
!     xknot(i) <= x < xknot(i+1)
!

    leftx = 0
    do ix = 1,nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(6,*) "subroutine dbsder:"
          write(6,*) "xknot(ix) <= xknot(ix+1) required."
          stop
       endif
       if ((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
    end do

    if (leftx .eq. 0) then
       write(6,*) "subroutine dbsder:"
       write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
       write(6,*) "xknot(1)     = ", xknot(1)
       write(6,*) "xknot(nx+kx) = ", xknot(nx+kx)
       write(6,*) "         x   = ", x
       stop
    endif

    if (iderx .eq. 0) then

       do ik = 1,kx-1
          work(ik) = bcoef(leftx+ik-kx)
          dl(ik)   = x - xknot(leftx+ik-kx)
          dr(ik)   = xknot(leftx+ik) - x
       end do

       work(kx)  = bcoef(leftx)
       dl(kx)    = x - xknot(leftx)

       do ik = 1,kx-1
          save2 = work(ik)
          do il = ik+1,kx
             save1 = work(il)
             work(il) = (dl(il) * work(il) + dr(il-ik) * save2)               &
                  &              / (dl(il) + dr(il - ik))
             save2 = save1
          end do
       end do

       dbsder = work(kx)

    elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then

       bsp(1) = 1.0_dabl
       do ik = 1,kx-iderx-1
          dr(ik) = xknot(leftx+ik) - x
          dl(ik) = x - xknot(leftx+1-ik)
          save   = bsp(1)
          bsp(1) = 0.0_dabl
          do il = 1, ik
             y         = save / (dr(il) + dl(ik+1-il))
             bsp(il)   = bsp(il) + dr(il) * y
             save      = bsp(il+1)
             bsp(il+1) = dl(ik+1-il) * y
          end do
       end do

       do ik = 1, kx
          work(ik) = bcoef(leftx+ik-kx)
          dr(ik)   = xknot(leftx+ik) - x
          dl(ik)   = x - xknot(leftx+ik-kx)
       end do

       do ik = 1, iderx
          dik   = dble(kx - ik)
          save2 = work(ik)
          do il = ik+1, kx
             save1    = work(il)
             work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik))
             save2    = save1
          end do
       end do

       sum = 0.0_dabl

       do ix = 1, kx-iderx
          sum = sum + bsp(ix) * work(iderx+ix)
       end do

       dbsder = sum

    else
       dbsder = 0.0_dabl
    endif
deallocate(work,dl,dr,bsp)

  end function dbsder


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbs1gd(iderx,nxvec,xvec,kx,xknot,nx,bcoef,val)

!
!  Evaluates the derivative of a spline on a grid, given its B-spline
!  representation.
!
!   iderx  - order of the derivative to be evaluated.  (input)
!            in particular, iderx = 0 returns the value of the
!            spline.
!   nxvec  - length of vector xvec.  (input)
!   xvec   - array of length nxvec containing the points at which the
!            spline is to be evaluated.  (input)
!            xvec should be strictly increasing.
!   kx     - order of the spline.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   nx     - number of B-spline coefficients.  (input)
!   bcoef  - array of length nx containing the B-spline
!            coefficients.  (input)
!   val    - array of length nxvec containing the values of the
!            iderx-th derivative of the spline at the points in
!            xvec.  (output)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                           :: iderx, nxvec, kx, nx
    real(kind=dabl), dimension(nxvec), intent(in)  :: xvec
    real(kind=dabl), dimension(nx), intent(in)     :: bcoef
    real(kind=dabl), dimension(nx+kx), intent(in)  :: xknot
    real(kind=dabl), dimension(nxvec), intent(out) :: val

    integer(kind=i4)                             :: i, il, ik, ix
    integer(kind=i4), dimension(:),allocatable           :: leftx
    real(kind=dabl)                      :: dik
    real(kind=dabl), dimension(:,:),allocatable :: dl, dr, biatx, work
    real(kind=dabl), dimension(:),allocatable    :: save1, save2, term

    logical :: same, next
allocate(leftx(nxvec),save1(nxvec), save2(nxvec), term(nxvec),dl(nxvec,kx), dr(nxvec,kx), biatx(nxvec,kx), work(nxvec,kx))

    leftx(1) = 0

    call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

    do ix = 2, nxvec
       leftx(ix) = leftx(ix-1)
       same = (xknot(leftx(ix)) .le. xvec(ix))                                &
            &        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
       if(.not. same ) then
          leftx(ix) = leftx(ix) + 1
          next      = (xknot(leftx(ix)) .le. xvec(ix))                        &
               &           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
          if (.not. next)                                                     &
               &           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
       endif
    end do

    do ix = 1, nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(6,*) "subroutine dbs1gd:"
          write(6,*) "xknot(ix) <= xknot(ix+1) required."
          write(6,*) ix, xknot(ix), xknot(ix+1)
          write(6,*)
          write(6,*) xknot
          stop
       endif
    end do

    do ix = 1, nxvec
       if ((xvec(ix).lt.xknot(1)).or.(xvec(ix).gt.xknot(nx+kx))) then
          write(6,*) "subroutine dbs1gd:"
          write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
          write(6,*) "x = ", xvec(ix)
          stop
       endif
    end do

    if (iderx .eq. 0) then

       do ix = 1,nxvec
          biatx(ix,1) = 1._dabl
          val(ix)     = 0._dabl
       end do

       do ik = 1, kx-1
          do ix = 1, nxvec
             dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
             dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
             save1(ix) = 0._dabl
          end do

          do il = 1, ik
             do ix = 1,nxvec
                term(ix)     = biatx(ix,il)                                   &
                     &                 / (dr(ix,il) + dl(ix,ik+1-il))
                biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                save1(ix)    = dl(ix,ik+1-il) * term(ix)
             end do
          end do

          do ix = 1, nxvec
             biatx(ix,ik+1) = save1(ix)
          end do
       end do

       do ik = 1, kx
          do ix = 1, nxvec
             val(ix) = val(ix) + biatx(ix,ik) * bcoef(leftx(ix)-kx+ik)
          end do
       end do

    elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then

       do ix = 1, nxvec
          biatx(ix,1) = 1._dabl
          val(ix)     = 0._dabl
       end do

       do ik = 1, kx-iderx-1
          do ix = 1, nxvec
             dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix)
             dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+1-ik)
             save1(ix)    = biatx(ix,1)
             biatx(ix,1) = 0.0_dabl
             do il = 1, ik
                term(ix)       = save1(ix)                                    &
                     &                 / (dr(ix,il) + dl(ix,ik+1-il))
                biatx(ix,il)   = biatx(ix,il) + dr(ix,il) * term(ix)
                save1(ix)      = biatx(ix,il+1)
                biatx(ix,il+1) = dl(ix,ik+1-il) * term(ix)
             end do
          end do
       end do

       do ik = 1, kx
          do ix = 1, nxvec
             work(ix,ik) = bcoef(leftx(ix)+ik-kx)
             dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix)
             dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+ik-kx)
          end do
       end do

       do ik = 1, iderx
          dik   = dble(kx - ik)
          do ix = 1, nxvec
             save2(ix) = work(ix,ik)
             do il = ik+1, kx
                save1(ix)   = work(ix,il)
                work(ix,il) = dik * (work(ix,il) - save2(ix))                 &
                     &                 /(dl(ix,il) + dr(ix,il-ik))
                save2(ix)   = save1(ix)
             end do
          end do
       end do

       do i = 1, kx-iderx
          do ix = 1, nxvec
             val(ix) = val(ix) + biatx(ix,i) * work(ix,iderx+i)
          end do
       end do

    else

       do ix = 1, nxvec
          val(ix) = 0.0_dabl
       end do

    endif
deallocate(leftx,save1, save2, term,dl, dr, biatx, work)

  end subroutine dbs1gd


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbsdca(iderx,x,kx,xknot,nx,bcoef,leftx)

!
! This routine is equivalent to the routine dbsder, but it does not
! check the parameters!!!
!
! Evaluates the derivative of a spline, given its B-spline representation.
!
!
!   iderx  - order of the derivative to be evaluated.  (input)
!            in particular, iderx = 0 returns the value of the
!            spline.
!   x      - point at which the spline is to be evaluated.  (input)
!   kx     - order of the spline.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   nx     - number of B-spline coefficients.  (input)
!   bcoef  - array of length nx containing the B-spline
!            coefficients.  (input)
!   leftx  - number of the intervall of xknot that includes x
!   dbsdca - value of the ideriv-th derivative of the spline at x.
!            (output)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                          :: iderx, kx, nx
    real(kind=dabl)                               :: dbsdca
    real(kind=dabl), intent(in)                   :: x
    real(kind=dabl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dabl), dimension(nx), intent(in)    :: bcoef

    integer(kind=i4)                       :: i, ik, il, leftx
    real(kind=dabl)                :: save, save1, save2, y, sum, dik
    real(kind=dabl), dimension(:),allocatable :: work, dl, dr,bsp
allocate(work(kx), dl(kx), dr(kx),bsp(kx))

    if (iderx .eq. 0) then

       do ik = 1, kx-1
          work(ik) = bcoef(leftx+ik-kx)
          dl(ik)   = x - xknot(leftx+ik-kx)
          dr(ik)   = xknot(leftx+ik) - x
       end do

       work(kx)  = bcoef(leftx)
       dl(kx)    = x - xknot(leftx)

       do ik = 1, kx-1
          save2 = work(ik)
          do il = ik+1, kx
             save1 = work(il)
             work(il) = (dl(il) * work(il) + dr(il-ik) * save2)               &
                  &              / (dl(il) + dr(il - ik))
             save2 = save1
          end do
       end do

       dbsdca = work(kx)

    elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then
       bsp(1) = 1.0_dabl
       do ik = 1,kx-iderx-1
          dr(ik) = xknot(leftx+ik) - x
          dl(ik) = x - xknot(leftx+1-ik)
          save   = bsp(1)
          bsp(1) = 0.0_dabl
          do il = 1, ik
             y         = save / (dr(il) + dl(ik+1-il))
             bsp(il)   = bsp(il) + dr(il) * y
             save      = bsp(il+1)
             bsp(il+1) = dl(ik+1-il) * y
          end do
       end do

       do ik = 1, kx
          work(ik) = bcoef(leftx+ik-kx)
          dr(ik)   = xknot(leftx+ik) - x
          dl(ik)   = x - xknot(leftx+ik-kx)
       end do

       do ik = 1, iderx
          dik   = dble(kx - ik)
          save2 = work(ik)
          do il = ik+1, kx
             save1    = work(il)
             work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik))
             save2    = save1
          end do
       end do

       sum = 0.0_dabl

       do i = 1, kx-iderx
          sum = sum + bsp(i) * work(iderx+i)
       end do

       dbsdca = sum

    else
       dbsdca = 0.0_dabl
    endif
deallocate(work, dl, dr,bsp)

  end function dbsdca


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbs2in(nx,xvec,ny,yvec,xydata,ldf,kx,ky,xknot,yknot,bcoef)

!
!  Computes a two-dimensional tensor-product spline interpolant,
!  returning the tensor-product B-spline coefficients.
!
!    nx     - number of data points in the x-direction.  (input)
!    xvec   - array of length nx containing the data points in
!             the x-direction.  (input)
!             xdata must be strictly increasing.
!    ny     - number of data points in the y-direction.  (input)
!    yvec   - array of length ny containing the data points in
!             the y-direction.  (input)
!             ydata must be strictly increasing.
!    xydata - array of size nx by nydata containing the values to
!             be interpolated.  (input)
!             fdata(i,j) is the value at (xdata(i),ydata(j)).
!    ldf    - the leading dimension of fdata exactly as specified in
!             the dimension statement of the calling program.
!             (input)
!    kx     - order of the spline in the x-direction.  (input)
!             kxord must be less than or equal to nxdata.
!    ky     - order of the spline in the y-direction.  (input)
!             kyord must be less than or equal to nydata.
!    xknot  - array of length nx+kx containing the knot
!             sequence in the x-direction.  (input)
!             xknot must be nondecreasing.
!    yknot  - array of length ny+ky containing the knot
!             sequence in the y-direction.  (input)
!             yknot must be nondecreasing.
!    bcoef  - array of length nx*ny containing the
!             tensor-product B-spline coefficients.  (output)
!             bscoef is treated internally as a matrix of size nxdata
!             by nydata.
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                           :: nx, ny, kx, ky, ldf

    real(kind=dabl), dimension(nx), intent(in)     :: xvec
    real(kind=dabl), dimension(ny), intent(in)     :: yvec
    real(kind=dabl), dimension(nx+kx), intent(in)  :: xknot
    real(kind=dabl), dimension(ny+ky), intent(in)  :: yknot
    real(kind=dabl), dimension(ldf,*), intent(in)  :: xydata
    real(kind=dabl), dimension(nx,ny), intent(out) :: bcoef

    real(kind=dabl), dimension(:,:),allocatable        :: work1
    real(kind=dabl), dimension(:),allocatable                   :: work2
    real(kind=dabl), dimension(:),allocatable :: work3
allocate(work1(max(nx,ny),max(nx,ny)),work2(max(nx,ny)),work3(max((2*kx-1)*nx,(2*ky-1)*ny)))

    call spli2d(xvec,ldf,xydata,xknot,nx,kx,ny,work2,work3,work1)
    call spli2d(yvec,ny, work1, yknot,ny,ky,nx,work2,work3,bcoef)
deallocate(work1,work2,work3)

  end subroutine dbs2in


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine spli2d(xyvec,ld,xydata,xyknot,n,k,m,work2,work3,bcoef)

    use numeric

    implicit none


    integer(kind=i4), intent(in)                         :: ld, n, k, m
    real(kind=dabl), dimension(n), intent(in)    :: xyvec
    real(kind=dabl), dimension(n+k), intent(in)  :: xyknot
    real(kind=dabl), dimension(ld,m), intent(in) :: xydata
    real(kind=dabl), dimension(m,n), intent(out) :: bcoef

    real(kind=dabl), dimension(n), intent(out)         :: work2
    real(kind=dabl), dimension((2*k-1)*n), intent(out) :: work3


    integer(kind=i4)        :: np1, km1, kpkm2, left, lenq, i, iflag, ilp1mx, j, jj
    real(kind=dabl) :: xyveci

    np1   = n + 1
    km1   = k - 1
    kpkm2 = 2 * km1
    left  = k
    lenq  = n * (k + km1)

    do i = 1,lenq
       work3(i) = 0.0_dabl
    end do

    do i = 1, n
       xyveci  = xyvec(i)
       ilp1mx = min0(i+k,np1)
       left   = max0(left,i)
       if (xyveci .lt. xyknot(left)) go to 998
30     if (xyveci .lt. xyknot(left+1)) go to 40
       left = left + 1
       if (left .lt. ilp1mx) go to 30
       left = left - 1
       if (xyveci .gt. xyknot(left+1)) go to 998
40     call bsplvb(xyknot,n+k,k,1,xyveci,left,work2)
       jj = i - left + 1 + (left - k) * (k + km1)
       do j = 1, k
          jj        = jj + kpkm2
          work3(jj) = work2(j)
       end do
    end do

    call banfac(work3,k+km1,n,km1,km1,iflag )

    if (iflag .ne. 1) then
       write(6,*) "subroutine dbs2in: error"
       write(6,*) "no solution of linear equation system !!!"
       stop
    end if

    do j = 1, m
       do i = 1, n
          work2(i) = xydata(i,j)
       end do

       call banslv(work3,k+km1,n,km1,km1,work2)

       do i = 1, n
          bcoef(j,i) = work2(i)
       end do
    end do

    return

998 write(6,*) "subroutine db2in:"
    write(6,*) "i with knot(i) <= x/y < knot(i+1) required."
    write(6,*) "knot(1)   = ", xyknot(1)
    write(6,*) "knot(n+k) = ", xyknot(n+k)
    write(6,*) "      x/y = ", xyveci

    stop

  end subroutine spli2d


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,bcoef)

!
!  evaluates a two-dimensional tensor-product spline, given its
!  tensor-product B-spline representation.    use numeric
!
!   x      - x-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   y      - y-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   bcoef  - array of length nx*ny containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny.
!   dbs2vl - value of the spline at (x,y).  (output)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                          :: nx, ny, kx, ky
    real(kind=dabl), intent(in)                   :: x, y
    real(kind=dabl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dabl), dimension(ny+ky), intent(in) :: yknot
    real(kind=dabl), dimension(nx,ny), intent(in) :: bcoef
    real(kind=dabl)                               :: dbs2vl

    integer(kind=i4)                       :: ix, iy, iky, leftx, lefty
    real(kind=dabl), dimension(:),allocatable :: work
allocate(work(ky))
!
!     check if knot(i) <= knot(i+1) and calculation of i so that
!     knot(i) <= x < knot(i+1)
!

    leftx = 0

    do ix = 1, nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(6,*) "subroutine dbs2vl:"
          write(6,*) "xknot(ix) <= xknot(ix+1) required."
          write(6,*) ix, xknot(ix), xknot(ix+1)
          write(6,*)
          write(6,*) xknot
          stop
       endif
       if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
    end do

    if(leftx .eq. 0) then
       write(6,*) "subroutine dbs2vl:"
       write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
       write(6,*) "x = ", x
       write(6,*)
       write(6,*) xknot
       stop
    endif

    lefty = 0

    do iy = 1, ny+ky-1
       if (yknot(iy) .gt. yknot(iy+1)) then
          write(6,*) "subroutine dbs2vl:"
          write(6,*) "yknot(iy) <= yknot(iy+1) required."
          write(6,*) iy, yknot(iy), yknot(iy+1)
          stop
       endif
       if((yknot(iy) .le. y) .and. (y .lt. yknot(iy+1))) lefty = iy
    end do

    if(lefty .eq. 0) then
       write(6,*) "subroutine dbs2vl:"
       write(6,*) "iy with yknot(iy) <= y < yknot(iy+1) required."
       write(6,*) "yknot(iy)   = ", yknot(iy)
       write(6,*) "  y         = ", y
       write(6,*) "yknot(iy+1) = ", yknot(iy+1)
       stop
    endif

    do iky = 1, ky
       work(iky) = dbsdca(0,x,kx,xknot,nx,bcoef(1,lefty-ky+iky),leftx)
    end do

    dbs2vl = dbsval(y,ky,yknot(lefty-ky+1),ky,work)
deallocate(work)
  end function dbs2vl


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,bcoef)

!
!  Evaluates the derivative of a two-dimensional tensor-product spline,
!  given its tensor-product B-spline representation.
!
!   iderx  - order of the derivative in the x-direction.  (input)
!   idery  - order of the derivative in the y-direction.  (input)
!   x      - x-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   y      - y-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   bcoef  - array of length nx*ny containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny.
!   dbs2dr  - value of the (iderx,idery) derivative of the spline at
!            (x,y).  (output)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                          :: iderx, idery
    integer(kind=i4), intent(in)                          :: kx, nx, ky, ny
    real(kind=dabl)                               :: dbs2dr
    real(kind=dabl), intent(in)                   :: x, y
    real(kind=dabl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dabl), dimension(ny+ky), intent(in) :: yknot
    real(kind=dabl), dimension(nx,ny), intent(in) :: bcoef

    integer(kind=i4)                       :: ix, iy, iky, nintx, ninty
    real(kind=dabl), dimension(:),allocatable :: work
allocate(work(ky))

!
!     check if knot(i) <= knot(i+1) and calculation of i so that
!     knot(i) <= x < knot(i+1)
!

    nintx = 0

    do ix = 1, nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(6,*) "subroutine dbs2dr:"
          write(6,*) "xknot(ix) <= xknot(ix+1) required."
          write(6,*) ix, xknot(ix), xknot(ix+1)
          stop
       endif
       if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) nintx = ix
    end do

    if(nintx .eq. 0) then
       write(6,*) "subroutine dbs2dr:"
       write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
       write(6,*) "x = ", x
       stop
    endif

    ninty = 0

    do iy = 1, ny+ky-1
       if (yknot(iy) .gt. yknot(iy+1)) then
          write(6,*) "subroutine dbs2dr:"
          write(6,*) "yknot(iy) <= yknot(iy+1) required."
          write(6,*) iy, yknot(iy), yknot(iy+1)
          stop
       endif
       if ((yknot(iy) .le. y) .and. (y .lt. yknot(iy+1))) ninty = iy
    end do

    if(ninty .eq. 0) then
       write(6,*) "subroutine dbs2dr:"
       write(6,*) "iy with yknot(iy) <= y < yknot(iy+1) required."
       write(6,*) "y = ", y
       stop
    endif

    do iky = 1, ky
       work(iky) =  dbsdca(iderx,x,kx,xknot,nx,bcoef(1,ninty-ky+iky),nintx)
    end do

    dbs2dr = dbsder(idery,y,ky,yknot(ninty-ky+1),ky,work)
deallocate(work)

  end function dbs2dr


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbs2gd(iderx,idery,nxvec,xvec,nyvec,yvec,kx,ky,xknot,yknot,      &
       & nx,ny,bcoef,val,ldf)

!
!  Evaluates the derivative of a two-dimensional tensor-product spline,
!  given its tensor-product B-spline representation on a grid.
!
!   iderx   - order of the derivative in the x-direction.  (input)
!   idery   - order of the derivative in the y-direction.  (input)
!   nxvec   - number of grid points in the x-direction.  (input)
!   xvec    - array of length nx containing the x-coordinates at
!             which the spline is to be evaluated.  (input)
!             the points in xvec should be strictly increasing.
!   nyvec   - number of grid points in the y-direction.  (input)
!   yvec    - array of length ny containing the y-coordinates at
!             which the spline is to be evaluated.  (input)
!             the points in yvec should be strictly increasing.
!   kx      - order of the spline in the x-direction.  (input)
!   ky      - order of the spline in the y-direction.  (input)
!   xknot   - array of length nx+kx containing the knot
!             sequence in the x-direction.  (input)
!             xknot must be nondecreasing.
!   yknot   - array of length ny+ky containing the knot
!             sequence in the y-direction.  (input)
!             yknot must be nondecreasing.
!   nx      - number of B-spline coefficients in the x-direction.
!             (input)
!   ny      - number of B-spline coefficients in the y-direction.
!             (input)
!   bcoef   - array of length nx*ny containing the
!             tensor-product B-spline coefficients.  (input)
!             bscoef is treated internally as a matrix of size nx
!             by ny.
!   val     - value of the (iderx,idery) derivative of the spline on
!             the nx by ny grid.  (output)
!             value(i,j) contains the derivative of the spline at the
!             point (xvec(i),yvec(j)).
!   ldf     - leading dimension of value exactly as specified in the
!             dimension statement of the calling program.  (input)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                           :: iderx, idery
    integer(kind=i4), intent(in)                           :: nxvec, nyvec
    integer(kind=i4), intent(in)                           :: kx, nx, ky, ny
    integer(kind=i4), intent(in)                           :: ldf

    real(kind=dabl), dimension(nxvec), intent(in)  :: xvec
    real(kind=dabl), dimension(nyvec), intent(in)  :: yvec
    real(kind=dabl), dimension(nx+kx), intent(in)  :: xknot
    real(kind=dabl), dimension(ny+ky), intent(in)  :: yknot
    real(kind=dabl), dimension(nx,ny), intent(in)  :: bcoef
    real(kind=dabl), dimension(ldf,*), intent(out) :: val

    integer(kind=i4)                                     :: i, ik, il, ix, iy, ikx, iky
    integer(kind=i4), dimension(:),allocatable                   :: leftx
    integer(kind=i4), dimension(:),allocatable                   :: lefty
    real(kind=dabl), dimension(:,:),allocatable         :: dl, dr
    real(kind=dabl), dimension(:),allocatable :: save1
    real(kind=dabl), dimension(:,:),allocatable         :: biatx
    real(kind=dabl), dimension(:,:),allocatable         :: biaty
    real(kind=dabl), dimension(:),allocatable :: term
    real(kind=dabl), dimension(:),allocatable               :: work

    logical :: same,next
allocate(leftx(nxvec),lefty(nyvec),dl(nxvec,kx),dr(nxvec,kx),save1(max(nxvec,nyvec)))
allocate(biatx(nxvec,kx),biaty(nyvec,ky),term(max(nxvec,nyvec)),work(ky))

    leftx(1) = 0

    call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

    do ix = 2, nxvec
       leftx(ix) = leftx(ix-1)
       same = (xknot(leftx(ix)) .le. xvec(ix))                                &
            &        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
       if(.not. same ) then
          leftx(ix) = leftx(ix) + 1
          next      = (xknot(leftx(ix)) .le. xvec(ix))                        &
               &           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
          if (.not. next)                                                     &
               &           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
       endif
    end do

    do i = 1, nx+kx-1
       if (xknot(i) .gt. xknot(i+1)) then
          write(6,*) "subroutine dbs2gd:"
          write(6,*) "xknot(i) <= xknot(i+1) required."
          write(6,*) i, xknot(i), xknot(i+1)
          write(6,*)
          write(6,*) xknot
          stop
       endif
    end do

    do i = 1, nxvec
       if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
          write(6,*) "subroutine dbs2gd:"
          write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
          write(6,*) "x = ", xvec(i)
          stop
       endif
    end do

    lefty(1) = 0

    call huntn(yknot,ny+ky,ky,yvec(1),lefty(1))

    do iy = 2, nyvec
       lefty(iy) = lefty(iy-1)
       same = (yknot(lefty(iy)) .le. yvec(iy))                                &
            &        .and. (yvec(iy) .le. yknot(lefty(iy)+1))
       if(.not. same ) then
          lefty(iy) = lefty(iy) + 1
          next      = (yknot(lefty(iy)) .le. yvec(iy))                        &
               &           .and. (yvec(iy) .le. yknot(lefty(iy)+1))
          if (.not. next) call huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy))
       endif
    end do

    do i = 1, ny+ky-1
       if (yknot(i) .gt. yknot(i+1)) then
          write(6,*) "subroutine dbs2gd:"
          write(6,*) "yknot(i) <= yknot(i+1) required."
          write(6,*) i, yknot(i), yknot(i+1)
          write(6,*)
          write(6,*) yknot
          stop
       endif
    end do

    do i = 1, nyvec
       if ((yvec(i).lt.yknot(1)).or.(yvec(i).gt.yknot(ny+ky))) then
          write(6,*) "subroutine dbs2gd:"
          write(6,*) "iy with yknot(iy) <= y < yknot(iy+1) required."
          write(6,*) "y = ", yvec(i)
          stop
       endif
    end do

    if ((iderx .eq. 0) .and. (idery .eq. 0)) then

       do ix = 1,nxvec
          biatx(ix,1) = 1._dabl
       end do

       do ik = 1, kx-1
          do ix = 1,nxvec
             dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
             dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
             save1(ix) = 0._dabl
          end do

          do il = 1,ik
             do ix = 1,nxvec
                term(ix)     = biatx(ix,il)                                   &
                     &                 / (dr(ix,il) + dl(ix,ik+1-il))
                biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                save1(ix)    = dl(ix,ik+1-il) * term(ix)
             end do
          end do

          do ix = 1, nxvec
             biatx(ix,ik+1) = save1(ix)
          end do
       end do

       do iy = 1, nyvec
          biaty(iy,1) = 1._dabl
       end do

       do ik = 1, ky-1
          do iy = 1, nyvec
             dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy)
             dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik)
             save1(iy) = 0._dabl
          end do

          do il = 1, ik
             do iy = 1,nyvec
                term(iy)     = biaty(iy,il)                                   &
                     &                 / (dr(iy,il) + dl(iy,ik+1-il))
                biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy)
                save1(iy)    = dl(iy,ik+1-il) * term(iy)
             end do
          end do

          do iy = 1, nyvec
             biaty(iy,ik+1) = save1(iy)
          end do
       end do

       do iy = 1, nyvec
          do ix = 1, nxvec
             val(ix,iy) = 0.0_dabl
          end do
       end do

       do iky = 1, ky
          do ikx = 1, kx
             do iy = 1, nyvec
                do ix = 1, nxvec
                   val(ix,iy) = val(ix,iy)                                    &
                        & + biatx(ix,ikx) * biaty(iy,iky)                     &
                        & * bcoef(leftx(ix)-kx+ikx,lefty(iy)-ky+iky)
                end do
             end do
          end do
       end do

    elseif (((iderx .ge. 1) .or. (idery .ge. 1))                              &
         &  .and. ( (iderx .lt. kx) .and. (idery .lt. ky))) then

       do iy = 1, nyvec
          do ix = 1, nxvec
             do iky = 1, ky
                work(iky) = dbsdca(iderx,xvec(ix),kx,xknot,nx,                &
                     &             bcoef(1,lefty(iy)-ky+iky),leftx(ix))
             end do
             val(ix,iy) = dbsder(idery,yvec(iy),ky,                           &
                  &              yknot(lefty(iy)-ky+1),ky,work)
          end do
       end do

    else

       do iy = 1, nyvec
          do ix = 1, nxvec
             val(ix,iy) = 0.0_dabl
          end do
       end do

    endif
deallocate(leftx,lefty,dl,dr,save1)
deallocate(biatx,biaty,term,work)

  end subroutine dbs2gd


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbs3in(nx,xvec,ny,yvec,nz,zvec,xyzdata,ldf,mdf,kx,ky,kz,         &
       & xknot,yknot,zknot,bcoef)

!
!  Computes a three-dimensional tensor-product spline interpolant,
!  returning the tensor-product B-spline coefficients.
!
!   nx      - number of data points in the x-direction.  (input)
!   xvec    - array of length nxdata containing the data points in
!             the x-direction.  (input)
!             xdata must be increasing.
!   ny      - number of data points in the y-direction.  (input)
!   yvec    - array of length nydata containing the data points in
!             the y-direction.  (input)
!             ydata must be increasing.
!   nz      - number of data points in the z-direction.  (input)
!   zvec    - array of length nzdata containing the data points in
!             the z-direction.  (input)
!             zdata must be increasing.
!   xyzdata - array of size nx by ny by nz containing the
!             values to be interpolated.  (input)
!             xyzdata(i,j,k) contains the value at
!             (xvec(i),yvec(j),zvec(k)).
!   ldf     - leading dimension of fdata exactly as specified in the
!             dimension statement of the calling program.  (input)
!   mdf     - middle dimension of fdata exactly as specified in the
!             dimension statement of the calling program.  (input)
!   kx      - order of the spline in the x-direction.  (input)
!             kxord must be less than or equal to nxdata.
!   ky      - order of the spline in the y-direction.  (input)
!             kyord must be less than or equal to nydata.
!   kz      - order of the spline in the z-direction.  (input)
!             kzord must be less than or equal to nzdata.
!   xknot   - array of length nx+kx containing the knot
!             sequence in the x-direction.  (input)
!             xknot must be nondecreasing.
!   yknot   - array of length ny+ky containing the knot
!             sequence in the y-direction.  (input)
!             yknot must be nondecreasing.
!   zknot   - array of length nz+kz containing the knot
!             sequence in the z-direction.  (input)
!             zknot must be nondecreasing.
!   bcoef   - array of length nx*ny*nz containing the
!             tensor-product B-spline coefficients.  (output)
!             bscoef is treated internally as a matrix of size nx
!             by ny by nz.
!

    use numeric

    implicit none

    integer(kind=i4), intent(in) :: nx, ny, nz, kx, ky, kz
    integer(kind=i4), intent(in) :: ldf, mdf

    real(kind=dabl), dimension(nx), intent(in)         :: xvec
    real(kind=dabl), dimension(ny), intent(in)         :: yvec
    real(kind=dabl), dimension(nz), intent(in)         :: zvec
    real(kind=dabl), dimension(nx+kx), intent(in)      :: xknot
    real(kind=dabl), dimension(ny+ky), intent(in)      :: yknot
    real(kind=dabl), dimension(nz+kz), intent(in)      :: zknot
    real(kind=dabl), dimension(ldf,mdf,nz), intent(in) :: xyzdata
    real(kind=dabl), dimension(nx,ny,nz), intent(out)  :: bcoef

    integer(kind=i4)                                :: iz
    real(kind=dabl), dimension(:,:,:),allocatable    :: work1
    real(kind=dabl), dimension(:),allocatable          :: work2
    real(kind=dabl), dimension(:),allocatable :: work3
allocate(work1(nx,ny,nz),work2(nz),work3((2*kz-1)*nz))

    call spli3d(zvec,ldf,mdf,xyzdata,zknot,nz,kz,nx,ny,work2,work3,work1,     &
         &     nx,ny,nz)

    do iz = 1, nz
       call dbs2in(nx,xvec,ny,yvec,work1(1,1,iz),nx,kx,ky,xknot,yknot,        &
            &        bcoef(1,1,iz))
    end do
deallocate(work1,work2,work3)

  end subroutine dbs3in


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine spli3d(xyzvec,ldf,mdf,xyzdata,xyzknot,n,k,m,l,work2,work3,       &
       & bcoef,nx,ny,nz)

    use numeric

    implicit none

    integer(kind=i4), intent(in)                               :: ldf, mdf, n, k, m, l
    integer(kind=i4), intent(in)                               :: nx, ny, nz
    real(kind=dabl), dimension(n), intent(in)          :: xyzvec
    real(kind=dabl), dimension(n+k), intent(in)        :: xyzknot
    real(kind=dabl), dimension(ldf,mdf,*), intent(in)  :: xyzdata
    real(kind=dabl), dimension(nx,ny,nz), intent(out)  :: bcoef
    real(kind=dabl), dimension(n), intent(out)         :: work2
    real(kind=dabl), dimension((2*k-1)*n), intent(out) :: work3

    integer(kind=i4)        :: np1, km1, kpkm2, left, lenq, i, ilp1mx, j, jj, iflag, in
    real(kind=dabl) :: xyzveci


    np1   = n + 1
    km1   = k - 1
    kpkm2 = 2 * km1
    left  = k
    lenq  = n * (k + km1)

    do i = 1, lenq
       work3(i) = 0._dabl
    end do

    do i = 1, n
       xyzveci = xyzvec(i)
       ilp1mx  = min0(i+k,np1)
       left    = max0(left,i)
       if (xyzveci .lt. xyzknot(left)) go to 998
30     if (xyzveci .lt. xyzknot(left+1)) go to 40
       left = left + 1
       if (left .lt. ilp1mx) go to 30
       left = left - 1
       if (xyzveci .gt. xyzknot(left+1)) go to 998
40     call bsplvb(xyzknot,n+k,k,1,xyzveci,left,work2)
       jj = i - left + 1 + (left - k) * (k + km1)
       do j = 1, k
          jj    = jj + kpkm2
          work3(jj) = work2(j)
       end do
    end do

    call banfac(work3,k+km1,n,km1,km1,iflag)

    if (iflag .ne. 1) then
       write(6,*) "subroutine dbs3in: error"
       write(6,*) "no solution of linear equation system !!!"
       stop
    end if

    do j = 1, l
       do i = 1, m
          do in = 1, n
             work2(in) = xyzdata(i,j,in)
          end do

          call banslv(work3,k+km1,n,km1,km1,work2)

          do in = 1, n
             bcoef(i,j,in) = work2(in)
          end do

       end do
    end do

    return

998 write(6,*) "subroutine db3in:"
    write(6,*) "i with knot(i) <= x/y/z < knot(i+1) required."
    write(6,*) "knot(1)   = ", xyzknot(1)
    write(6,*) "knot(n+k) = ", xyzknot(n+k)
    write(6,*) "    x/y/z = ", xyzveci

    stop

  end subroutine spli3d


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbs3vl(x,y,z,kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef)

!
!  Evaluates a three-dimensional tensor-product spline, given its
!  tensor-product B-spline representation.
!
!   x      - x-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   y      - y-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   z      - z-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   kz     - order of the spline in the z-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   zknot  - array of length nz+kz containing the knot
!            sequence in the z-direction.  (input)
!            zknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   nz     - number of B-spline coefficients in the z-direction.
!            (input)
!   bcoef  - array of length nx*ny*nz containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny by nz.
!   dbs3vl - value of the spline at (x,y,z).  (output)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                             :: nx, ny, nz, kx, ky, kz
    real(kind=dabl), intent(in)                      :: x, y, z
    real(kind=dabl), dimension(nx+kx), intent(in)    :: xknot
    real(kind=dabl), dimension(ny+ky), intent(in)    :: yknot
    real(kind=dabl), dimension(nz+kz), intent(in)    :: zknot
    real(kind=dabl), dimension(nx,ny,nz), intent(in) :: bcoef
    real(kind=dabl)                                  :: dbs3vl

    integer(kind=i4)                       :: iz, nintz
    real(kind=dabl), dimension(:),allocatable :: work
allocate(work(kz))
!
!     check if knot(i) <= knot(i+1) and calculation of i so that
!     knot(i) <= x < knot(i+1)
!

    nintz = 0

    do iz = 1, nz+kz-1
       if (zknot(iz) .gt. zknot(iz + 1)) then
          write(6,*) "subroutine dbs3vl:"
          write(6,*) "zknot(iz) <= zknot(iz+1) required."
          write(6,*) iz, zknot(iz), zknot(iz+1)
          stop
       endif
       if((zknot(iz) .le. z) .and. (z .lt. zknot(iz + 1))) nintz = iz
    end do

    if(nintz .eq. 0) then
       write(6,*) "subroutine dbs3vl:"
       write(6,*) "iz with zknot(iz) <= z < zknot(iz+1) required."
       write(6,*) "zknot(iz)   = ", zknot(iz)
       write(6,*) "  z         = ", z
       write(6,*) "zknot(iz+1) = ", zknot(iz+1)
       stop
    endif

    do iz = 1, kz
       work(iz) = dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,bcoef(1,1,nintz-kz+iz))
    end do

    dbs3vl = dbsval(z,kz,zknot(nintz-kz+1),kz,work)
deallocate(work)
  end function dbs3vl


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbs3dr(iderx,idery,iderz,x,y,z,kx,ky,kz,xknot,yknot,zknot,         &
       & nx,ny,nz,bcoef)

!
!  Evaluates the derivative of a three-dimensional tensor-product spline,
!  given its tensor-product B-spline representation.
!
!   iderx  - order of the x-derivative.  (input)
!   idery  - order of the y-derivative.  (input)
!   iderz  - order of the z-derivative.  (input)
!   x      - x-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   y      - y-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   z      - z-coordinate of the point at which the spline is to be
!            evaluated.  (input)
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   kz     - order of the spline in the z-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   zknot  - array of length nz+kz containing the knot
!            sequence in the z-direction.  (input)
!            zknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   nz     - number of B-spline coefficients in the z-direction.
!            (input)
!   bcoef  - array of length nx*ny*nz containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny by nz.
!   dbs3dr - value of the (iderx,idery,iderz) derivative of the
!            spline at (x,y,z).  (output)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                              :: iderx, idery, iderz
    integer(kind=i4), intent(in)                              :: nx, ny, nz, kx, ky, kz
    real(kind=dabl), intent(in)                       :: x, y, z
    real(kind=dabl), dimension(nx+kx), intent(in)     :: xknot
    real(kind=dabl), dimension(ny+ky), intent(in)     :: yknot
    real(kind=dabl), dimension(nz+kz), intent(in)     :: zknot
    real(kind=dabl), dimension(nx,ny,nz), intent(in)  :: bcoef
    real(kind=dabl)                                   :: dbs3dr

    integer(kind=i4)                       :: iz, nintz
    real(kind=dabl), dimension(:),allocatable :: work
allocate(work(kz))

!
!     check if knot(i) <= knot(i+1) and calculation of i so that
!     knot(i) <= x < knot(i+1)
!

    nintz = 0

    do iz = 1, nz+kz-1
       if (zknot(iz) .gt. zknot(iz + 1)) then
          write(6,*) "subroutine dbs3vl:"
          write(6,*) "zknot(iz) <= zknot(iz+1) required."
          write(6,*) iz, zknot(iz), zknot(iz+1)
          stop
       endif
       if((zknot(iz) .le. z) .and. (z .lt. zknot(iz + 1))) nintz = iz
    end do

    if(nintz .eq. 0) then
       write(6,*) "subroutine dbs3dr:"
       write(6,*) "iz with zknot(iz) <= z < zknot(iz+1) required."
       write(6,*) "zknot(iz)   = ", zknot(iz)
       write(6,*) "  z         = ", z
       write(6,*) "zknot(iz+1) = ", zknot(iz+1)
       stop
    endif

    do iz = 1, kz
       work(iz) = dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,             &
            &        bcoef(1,1,nintz-kz+iz))
    end do

    dbs3dr = dbsder(iderz,z,kz,zknot(nintz-kz+1),kz,work)
deallocate(work)

  end function dbs3dr


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbs3gd(iderx,idery,iderz,nxvec,xvec,nyvec,yvec,nzvec,zvec,       &
       & kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef,val,ldf,mdf)

!
!  Evaluates the derivative of a three-dimensional tensor-product spline,
!  given its tensor-product B-spline representation on a grid.
!
!   iderx  - order of the x-derivative.  (input)
!   idery  - order of the y-derivative.  (input)
!   iderz  - order of the z-derivative.  (input)
!   nx     - number of grid points in the x-direction.  (input)
!   xvec   - array of length nx containing the x-coordinates at
!            which the spline is to be evaluated.  (input)
!            the points in xvec should be strictly increasing.
!   ny     - number of grid points in the y-direction.  (input)
!   yvec   - array of length ny containing the y-coordinates at
!            which the spline is to be evaluated.  (input)
!            the points in yvec should be strictly increasing.
!   nz     - number of grid points in the z-direction.  (input)
!   zvec   - array of length nz containing the z-coordinates at
!            which the spline is to be evaluated.  (input)
!            the points in yvec should be strictly increasing.
!   kx     - order of the spline in the x-direction.  (input)
!   ky     - order of the spline in the y-direction.  (input)
!   kz     - order of the spline in the z-direction.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence in the x-direction.  (input)
!            xknot must be nondecreasing.
!   yknot  - array of length ny+ky containing the knot
!            sequence in the y-direction.  (input)
!            yknot must be nondecreasing.
!   zknot  - array of length nz+kz containing the knot
!            sequence in the z-direction.  (input)
!            zknot must be nondecreasing.
!   nx     - number of B-spline coefficients in the x-direction.
!            (input)
!   ny     - number of B-spline coefficients in the y-direction.
!            (input)
!   nz     - number of B-spline coefficients in the z-direction.
!            (input)
!   bcoef  - array of length nx*ny*nz containing the
!            tensor-product B-spline coefficients.  (input)
!            bscoef is treated internally as a matrix of size nx
!            by ny by nz.
!   val    - array of size nx by ny by nz containing the values of
!            the (iderx,idery,iderz) derivative of the spline on the
!            nx by ny by nz grid.  (output)
!            value(i,j,k) contains the derivative of the spline at
!            the point (xvec(i), yvec(j), zvec(k)).
!   ldf    - leading dimension of value exactly as specified in the
!            dimension statement of the calling program.  (input)
!   mdf    - middle dimension of value exactly as specified in the
!            dimension statement of the calling program.  (input)
!

    use numeric

    implicit none

    integer(kind=i4), intent(in)                               :: iderx, idery, iderz
    integer(kind=i4), intent(in)                               :: nxvec, nyvec, nzvec
    integer(kind=i4), intent(in)                               :: kx, nx, ky, ny, kz, nz
    integer(kind=i4), intent(in)                               :: ldf,mdf

    real(kind=dabl), dimension(nxvec), intent(in)      :: xvec
    real(kind=dabl), dimension(nyvec), intent(in)      :: yvec
    real(kind=dabl), dimension(nzvec), intent(in)      :: zvec
    real(kind=dabl), dimension(nx+kx), intent(in)      :: xknot
    real(kind=dabl), dimension(ny+ky), intent(in)      :: yknot
    real(kind=dabl), dimension(nz+kz), intent(in)      :: zknot
    real(kind=dabl), dimension(nx,ny,nz), intent(in)   :: bcoef
    real(kind=dabl), dimension(ldf,mdf,*), intent(out) :: val

    integer(kind=i4)                                           :: i, ik, il, ix, iy, iz
    integer(kind=i4)                                           :: ikx, iky, ikz
    integer(kind=i4), dimension(:),allocatable                         :: leftx
    integer(kind=i4), dimension(:),allocatable                         :: lefty
    integer(kind=i4), dimension(:),allocatable                         :: leftz
    real(kind=dabl), dimension(:,:),allocatable               :: biatx
    real(kind=dabl), dimension(:,:),allocatable               :: biaty
    real(kind=dabl), dimension(:,:),allocatable               :: biatz
    real(kind=dabl), dimension(:),allocatable :: term, save1

    real(kind=dabl), dimension(:,:),allocatable :: dl, dr

    logical :: same,next

allocate(leftx(nxvec),lefty(nyvec),leftz(nzvec),dl(max(nxvec,nyvec,nzvec), max(kx,ky,kz)),&
       & dr(max(nxvec,nyvec,nzvec), max(kx,ky,kz)),save1(max(nxvec,nyvec,nzvec)))
allocate(biatx(nxvec,kx),biaty(nyvec,ky),biatz(nzvec,kz),term(max(nxvec,nyvec,nzvec)))

    do i = 1, nx+kx-1
       if (xknot(i) .gt. xknot(i+1)) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "xknot(i) <= xknot(i+1) required."
          write(6,*) i, xknot(i), xknot(i+1)
          write(6,*)
          write(6,*) xknot
          stop
       endif
    end do

    do i = 1, nxvec
       if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
          write(6,*) "x = ", xvec(i)
          stop
       endif
    end do

    leftx(1) = 0

    call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

    do ix = 2, nxvec
       leftx(ix) = leftx(ix-1)
       same = (xknot(leftx(ix)) .le. xvec(ix))                                &
            &        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
       if(.not. same ) then
          leftx(ix) = leftx(ix) + 1
          next      = (xknot(leftx(ix)) .le. xvec(ix))                        &
               &           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
          if (.not. next) call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
       endif
    end do

    do i = 1, ny+ky-1
       if (yknot(i) .gt. yknot(i+1)) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "yknot(i) <= yknot(i+1) required."
          write(6,*) i, yknot(i), yknot(i+1)
          write(6,*)
          write(6,*) yknot
          stop
       endif
    end do

    do i = 1, nyvec
       if ((yvec(i).lt.yknot(1)).or.(yvec(i).gt.yknot(ny+ky))) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "iy with yknot(iy) <= y < yknot(iy+1) required."
          write(6,*) "y = ", yvec(i)
          stop
       endif
    end do

    lefty(1) = 0

    call huntn(yknot,ny+ky,ky,yvec(1),lefty(1))

    do iy = 2, nyvec
       lefty(iy) = lefty(iy-1)
       same = (yknot(lefty(iy)) .le. yvec(iy))                                &
            &        .and. (yvec(iy) .le. yknot(lefty(iy)+1))
       if(.not. same ) then
          lefty(iy) = lefty(iy) + 1
          next      = (yknot(lefty(iy)) .le. yvec(iy))                        &
               &           .and. (yvec(iy) .le. yknot(lefty(iy)+1))
          if (.not. next) call huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy))
       endif
    end do

    do i = 1,nz+kz-1
       if (zknot(i) .gt. zknot(i+1)) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "zknot(i) <= zknot(i+1) required."
          write(6,*) i, zknot(i), zknot(i+1)
          write(6,*)
          write(6,*) zknot
          stop
       endif
    end do

    do i = 1, nzvec
       if ((zvec(i).lt.zknot(1)).or.(zvec(i).gt.zknot(nz+kz))) then
          write(6,*) "subroutine dbs3gd:"
          write(6,*) "iz with zknot(iz) <= z < zknot(iz+1) required."
          write(6,*) "z = ", zvec(i)
          stop
       endif
    end do

    leftz(1) = 0

    call huntn(zknot,nz+kz,kz,zvec(1),leftz(1))

    do iz = 2, nzvec
       leftz(iz) = leftz(iz-1)
       same = (zknot(leftz(iz)) .le. zvec(iz))                                &
            &        .and. (zvec(iz) .le. zknot(leftz(iz)+1))
       if(.not. same ) then
          leftz(iz) = leftz(iz) + 1
          next      = (zknot(leftz(iz)) .le. zvec(iz))                        &
               &           .and. (zvec(iz) .le. zknot(leftz(iz)+1))
          if (.not. next) call huntn(zknot,nz+kz,kz,zvec(iz),leftz(iz))
       endif
    end do

    if ((iderx .eq. 0) .and. (idery .eq. 0) .and. (iderz .eq.0)) then

       do ix = 1, nxvec
          biatx(ix,1) = 1.0_dabl
       end do

       do ik = 1, kx-1
          do ix = 1, nxvec
             dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
             dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
             save1(ix) = 0._dabl
          end do

          do il = 1, ik
             do ix = 1, nxvec
                term(ix)     = biatx(ix,il) / (dr(ix,il) + dl(ix,ik+1-il))
                biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                save1(ix)    = dl(ix,ik+1-il) * term(ix)
             end do
          end do

          do ix = 1, nxvec
             biatx(ix,ik+1) = save1(ix)
          end do
       end do

       do iy = 1, nyvec
          biaty(iy,1) = 1.0_dabl
       end do

       do ik = 1, ky-1
          do iy = 1, nyvec
             dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy)
             dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik)
             save1(iy) = 0._dabl
          end do

          do il = 1,ik
             do iy = 1,nyvec
                term(iy)     = biaty(iy,il) / (dr(iy,il) + dl(iy,ik+1-il))
                biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy)
                save1(iy)    = dl(iy,ik+1-il) * term(iy)
             end do
          end do

          do iy = 1,nyvec
             biaty(iy,ik+1) = save1(iy)
          end do
       end do

       do iz = 1,nzvec
          biatz(iz,1) = 1.0_dabl
       end do

       do ik = 1, kz-1
          do iz = 1, nzvec
             dr(iz,ik) = zknot(leftz(iz)+ik) - zvec(iz)
             dl(iz,ik) = zvec(iz) - zknot(leftz(iz)+1-ik)
             save1(iz) = 0._dabl
          end do

          do il = 1, ik
             do iz = 1, nzvec
                term(iz)     = biatz(iz,il) / (dr(iz,il) + dl(iz,ik+1-il))
                biatz(iz,il) = save1(iz) + dr(iz,il) * term(iz)
                save1(iz)    = dl(iz,ik+1-il) * term(iz)
             end do
          end do

          do iz = 1, nzvec
             biatz(iz,ik+1) = save1(iz)
          end do
       end do

       do iz = 1,nzvec
          do iy = 1,nyvec
             do ix = 1,nxvec
                val(ix,iy,iz) = 0.0_dabl
             end do
          end do
       end do

       do ikz = 1, kz
          do iky = 1, ky
             do ikx = 1, kx
                do iz = 1, nzvec
                   do iy = 1, nyvec
                      do ix = 1, nxvec
                         val(ix,iy,iz) = val(ix,iy,iz)                        &
                              &  + biatx(ix,ikx) * biaty(iy,iky)              &
                              &  * biatz(iz,ikz)                              &
                              &  * bcoef(leftx(ix)-kx+ikx,                    &
                              &          lefty(iy)-ky+iky,leftz(iz)-kz+ikz)
                      end do
                   end do
                end do
             end do
          end do
       end do

    else

       do iz = 1, nzvec
          do iy = 1, nyvec
             do ix = 1, nxvec
                val(ix,iy,iz) = dbs3dr(iderx,idery,iderz,xvec(ix),            &
                     &  yvec(iy),zvec(iz),kx,ky,kz,xknot,yknot,               &
                     &  zknot,nx,ny,nz,bcoef)
             end do
          end do
       end do

    endif
deallocate(leftx,lefty,leftz,dl,&
       & dr,save1)
deallocate(biatx,biaty,biatz,term)

  end subroutine dbs3gd


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine bsplvb(t,n,jhigh,index,x,left,biatx)

    use numeric

    implicit none

    integer(kind=i4), intent(in) :: n, jhigh, index, left

    real(kind=dabl), intent(in)                    :: x
    real(kind=dabl), dimension(n), intent(in)      :: t
    real(kind=dabl), dimension(jhigh), intent(out) :: biatx

    integer(kind=i4)                          :: j = 1
    integer(kind=i4)                          :: i, jp1
    real(kind=dabl)                   :: saved, term
    real(kind=dabl), dimension(:),allocatable :: dl, dr
allocate(dl(jhigh),dr(jhigh))

    if (index .eq. 1) then
       j = 1
       biatx(1) = 1.0_dabl
       if (j .ge. jhigh) return
    end if

20  jp1 = j + 1

    dr(j) = t(left+j) - x
    dl(j) = x - t(left+1-j)
    saved = 0._dabl

    do i = 1, j
       term     = biatx(i) / (dr(i) + dl(jp1-i))
       biatx(i) = saved + dr(i) * term
       saved    = dl(jp1-i) * term
    end do

    biatx(jp1) = saved
    j          = jp1

    if (j .lt. jhigh) go to 20
deallocate(dl,dr)

  end subroutine bsplvb


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine banfac(w,nroww,nrow,nbandl,nbandu,iflag)

    use numeric

    implicit none

    integer(kind=i4), intent(in)                                  :: nroww,nrow
    integer(kind=i4), intent(in)                                  :: nbandl,nbandu
    integer(kind=i4), intent(out)                                 :: iflag
    real(kind=dabl), dimension(nroww,nrow), intent(inout) :: w

    real(kind=dabl) :: pivot, factor
    integer(kind=i4)        :: middle, nrowm1, jmax, kmax, ipk, midmk, i, j, k


    iflag  = 1
    middle = nbandu + 1
    nrowm1 = nrow - 1

    if (nrowm1 .lt. 0) goto 999
    if (nrowm1 .eq. 0) goto 900
    if (nrowm1 .gt. 0) goto 10

10  if (nbandl .gt. 0) go to 30

    do i = 1, nrowm1
       if (w(middle,i) .eq. 0._dabl) go to 999
    end do

    go to 900

30  if (nbandu .gt. 0) go to 60

    do i = 1, nrowm1
       pivot = w(middle,i)
       if(pivot .eq. 0._dabl) go to 999
       jmax = min0(nbandl, nrow - i)
       do j = 1, jmax
          w(middle+j,i) = w(middle+j,i) / pivot
       end do
    end do

    return

60  do i = 1, nrowm1
       pivot = w(middle,i)
       if (pivot .eq. 0._dabl) go to 999
       jmax = min0(nbandl,nrow - i)
       do j = 1,jmax
          w(middle+j,i) = w(middle+j,i) / pivot
       end do

       kmax = min0(nbandu,nrow - i)

       do k = 1, kmax
          ipk    = i + k
          midmk  = middle - k
          factor = w(midmk,ipk)
          do j = 1, jmax
             w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)                  &
                  &              * factor
          end do
       end do
    end do

900 if (w(middle,nrow) .ne. 0._dabl) return
999 iflag = 2

  end subroutine banfac


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine banslv(w,nroww,nrow,nbandl,nbandu,b)

    use numeric

    implicit none

    integer(kind=i4), intent(in)                               :: nroww,nrow
    integer(kind=i4), intent(in)                               :: nbandl,nbandu
    real(kind=dabl), dimension(nroww,nrow), intent(in) :: w
    real(kind=dabl), dimension(nrow), intent(inout)    :: b

    integer(kind=i4) :: middle, nrowm1, jmax, i, j

    middle = nbandu + 1
    if (nrow .eq. 1) goto 99
    nrowm1 = nrow - 1
    if (nbandl .eq. 0) goto 30

    do i = 1, nrowm1
       jmax = min0(nbandl, nrow - i)
       do j = 1, jmax
          b(i+j) = b(i+j) - b(i) * w(middle+j,i)
       end do
    end do

30  if (nbandu .gt. 0)  goto 50

    do i = 1, nrow
       b(i) = b(i) / w(1,i)
    end do

    return

50  do i = nrow, 2, -1
       b(i) = b(i)/w(middle,i)
       jmax = min0(nbandu,i-1)
       do j = 1, jmax
          b(i-j) = b(i-j) - b(i) * w(middle-j,i)
       end do
    end do

99  b(1) = b(1) / w(middle,1)

  end subroutine banslv


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine huntn(xx,n,kord,x,jlo)

    use numeric

    implicit none

    integer(kind=i4), intent(in)                      :: n, kord
    real(kind=dabl), intent(in)               :: x
    real(kind=dabl), dimension(n), intent(in) :: xx

    integer(kind=i4), intent(out)                     :: jlo

    integer(kind=i4) :: max, null, jhi, jm, inc

!
!     works only for B-Splines (order n)
!

    max  = n - kord
    null = kord

    if (jlo.le.null.or.jlo.gt.max) then
       jlo = null
       jhi = max+1
       goto 30
    endif

    inc = 1

    if (x .ge. xx(jlo)) then
10     jhi = jlo + inc
       if (jhi .gt. max) then
          jhi = max + 1
       else if (x .ge. xx(jhi)) then
          jlo = jhi
          inc = inc + inc
          goto 10
       endif
    else
       jhi = jlo
20     jlo = jhi - inc
       if (jlo .le. null) then
          jlo = null
       else if (x .lt. xx(jlo)) then
          jhi = jlo
          inc = inc + inc
          goto 20
       endif
    endif

30  if (jhi-jlo.eq.1) return

    jm = (jhi + jlo) / 2
    if (x .gt. xx(jm)) then
       jlo = jm
    else
       jhi = jm
    endif

    goto 30

  end subroutine huntn


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end module bspline
