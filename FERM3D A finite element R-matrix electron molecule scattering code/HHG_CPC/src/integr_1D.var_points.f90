!!!*************************************************************
! 文件/File: integr_1D.var_points.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: integr_1D.var_points.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************


SUBROUTINE qgaus(func,a,b,ss,x,w,ngp) 
  use nrtype
  REAL(kind=dbl) a,b,ss,func 
  integer(kind=i4b),intent(in) :: ngp
  EXTERNAL func 
  !Returns as ss the integral of the function func between a and b, by ten-point Gauss- Legendre
  !integration: the function is evaluated exactly ten times at interior points in
  !the range of integration. 
  REAL(kind=dbl) dx,xm,xr
  REAL(kind=dbl),intent(in) :: w(ngp),x(ngp)
  INTEGER j
  ! SAVE w,x 
  xm=0.5d0*(b+a)
  xr=0.5d0*(b-a) 
  ss=0
  do  j=1,ngp
     dx=xr*x(j)
     ss=ss+w(j)*(func(xm+dx)+func(xm-dx)) 
  enddo
  ss=xr*ss/2.d0 
  return 
END SUBROUTINE qgaus
SUBROUTINE  gauleg(ngp, xabsc, weig)
  use nrtype
  implicit none
  INTEGER  i, j, m
  REAL(dbl)  p1, p2, p3, pp, z, z1,m_Pi
  INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
  REAL(dbl), INTENT(OUT) :: xabsc(ngp), weig(ngp)
  real(kind=dbl),parameter :: eps=1.d-15
  m_Pi=Pi
  m = (ngp + 1) / 2
  !* Roots are symmetric in the interval - so only need to find half of them*/
  do i = 1, m                          ! Loop over the desired roots*/
     z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
     !*   Starting with the above approximation to the ith root,
     !*          we enter the main loop of refinement by NEWTON'S method   */
100  p1 = 1.0d0
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
     xabsc(i) =  - z                         ! Roots will be bewteen -1.0 & 1.0 */
     xabsc(ngp+1-i) =  + z                   ! and symmetric about theorigin */
     weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its*/
     weig(ngp+1-i) = weig(i)               ! symmetric counterpart*/
  end do     ! i loop
End subroutine gauleg
