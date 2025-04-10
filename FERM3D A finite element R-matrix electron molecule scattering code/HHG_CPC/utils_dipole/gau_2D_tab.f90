!!!*************************************************************
! 文件/File: gau_2D_tab.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: gau_2D_tab.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


module gauss2dint
  implicit none
  !  INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
  INTEGER,parameter,private :: dbp = kind(1.0d0)
  !   PRIVATE
  REAL (dbp),private :: newv
  REAL(dbp),private  :: EPS, M_PI
  PARAMETER (EPS=3.0d-15)              !EPS is the relative precision
  !  PARAMETER (M_PI=3.141592654d0)      ! Pi value
  PARAMETER (M_PI=3.1415926535897932384626433832d0)
  integer countx,county,counter

  !------------------------
  ! 2 dimensional integration subroutines
  !------------------------

contains

  RECURSIVE function qgss2d(origfn, xx1, xx2, yf1, yf2, ngp) RESULT(inth)
    implicit none                               ! returns integral in inth
    REAL(dbp)  inth, xx1, xx2, yf1, yf2, origfn
    REAL(dbp)  xm, xl, xtmp, ytmp
    INTEGER j
    INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
    REAL(dbp) :: xabsc(ngp), weig(ngp)
    external origfn

    !   interface
    !real(dbp)       function origfn(xp,yp,zp)
    !RESULT(vfun2d)     ! Original Function's Interface

    !                implicit none
    !         INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
    !                REAL(dbp)  vfun2d, xp, yp, zp
    !                end function origfn
    !      function yf1(x) RESULT(vy1)
    !        implicit none
    !        INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
    !        REAL(dbp)  vy1, x
    !      end  function yf1
    !      function yf2(x) RESULT(vy2)
    !        implicit none
    !        INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
    !        REAL(dbp)  vy2, x
    !     end  function yf2
    !      function zf1(x,y) RESULT(vz1)
    !        implicit none
    !        INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
    !        REAL(dbp)  vz1, x, y
    !      end  function zf1
    !      function zf2(x,y) RESULT(vz2)
    !        implicit none
    !        INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
    !        REAL(dbp)  vz2, x, y
    !      end  function zf2
    !   end interface
    call gauleg(ngp, xabsc, weig)
    inth = 0.0d0
    xm = 0.5d0 * (xx2 + xx1)
    xl = 0.5d0 * (xx2 - xx1)
    do j = 1, ngp
       xtmp = xm + xl*xabsc(j)                ! Gauss-Legendre Abcissas
       countx=countx+1
       county=0
       inth = inth + weig(j) * qgssgy()
    END do
    inth = inth * xl;    !Scale the answer to the range of integration  */
CONTAINS
    RECURSIVE function qgssgy() RESULT(intg)
      implicit none                                ! returns integral in intg
      REAL(dbp)  intg
      REAL(dbp)  ym, yl                                 ! all undeclared variables are
      INTEGER j                                    !   COOMON with HOST
      intg = 0.0d0
      ym = 0.5d0 * (yf2 + yf1)
      yl = 0.5d0 * (yf2 - yf1)
      do j = 1, ngp
         ytmp = ym + yl*xabsc(j)                ! Gauss-Legendre Abcissas
         county=county+1
         intg = intg + weig(j) * origfn(xtmp,ytmp)
      END do
      intg = intg * yl;    !Scale the answer to the range of integration  */
    END function qgssgy
  END FUNCTION  qgss2d


  SUBROUTINE  gauleg(ngp, xabsc, weig)
    implicit none
    INTEGER  i, j, m
    REAL(dbp)  p1, p2, p3, pp, z, z1
    INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
    REAL(dbp), INTENT(OUT) :: xabsc(ngp), weig(ngp)
    m = (ngp + 1) / 2
    !* Roots are symmetric in the interval - so only need to find half of them  */
    do i = 1, m                          ! Loop over the desired roots*/
       z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
       !*   Starting with the above approximation to the ith root,
       !*          we enter the main loop of refinement by NEWTON'S method   */
100    p1 = 1.0d0
       p2 = 0.0d0
       !*  Loop up the recurrence relation to get the Legendre
       !*  polynomial evaluated at z                 */
       do j = 1, ngp
          p3 = p2
          p2 = p1
          p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
       enddo
       !* p1 is now the desired Legendre polynomial. We next compute pp,
       !* its derivative, by a standard relation involving also p2, the
       !* polynomial of one lower order.      */
       pp = ngp*(z*p1-p2)/(z*z-1.0d0)
       z1 = z
       z = z1 - p1/pp             ! Newton's Method  */
       if (dabs(z-z1) .gt. EPS) GOTO  100
       xabsc(i) =  - z                         ! Roots will be bewteen -1.0 &1.0 */
       xabsc(ngp+1-i) =  + z                   ! and symmetric about theorigin  */
       weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its*/
       weig(ngp+1-i) = weig(i)               ! symmetric counterpart*/
    end do     ! i loop
  End subroutine gauleg
end module gauss2dint
