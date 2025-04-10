!!!*************************************************************
! 文件/File: vfitt.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vfitt.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

subroutine  vfitt(num_vib_levels)    
  ! this subroutine handles vlambda fitting at long and short range.
  !   we consider three different regions:
  !     (set nlast=npun-nfin)
  !     - a/r + b*r  short range region (i)
  !                 rbegin     ---> rr(2)
  !     - spline-interpolated medium range region (ii)
  !                 rr(2)      ---> rr(nlast-1)
  !     - inverse-power [c*r^(-n)+d*r^(-m)] extrapolated long range region (iii)
  !                 rr(nlast-1)---> rend
  !     (matching of the first derivatives at the boundary is ensured)
  use precision
  use shared, only: debug
  !integer debug 
  
  use potential, n_lambda=>num_lambda_terms, np=>num_pot_points 
  !integer :: num_lambda_terms
  !integer :: num_pot_points 
  !integer :: points_to_discard 
  !integer :: lambda_to_zero 
 
  !These are already allocated in main
  !real dimension(:,:,:,:) :: v_lambda, spline_coeff
  !real dimension(:)       :: r
  !real dimension(4,:,:,:) :: extrapolation_coeff

  use long_range_mod, n=>n_exp, m=>m_exp
  !integer, allocatable, dimension(:) :: n_exp, m_exp

  implicit none
  
  integer, intent(IN) :: num_vib_levels
  
  integer :: ir, l, n1, n2, iwrite, unit  
  integer :: n_last !last point to consider
  real(kind=wp), dimension(np) :: v_temp, v2_temp
  real(kind=wp) :: r_start, r_last, r_modal
  real(kind=wp) :: v_last, v_start
  real(kind=wp) :: a,b,c,d 
  real(kind=wp) :: short_range_derivative, long_range_derivative

  r_start=r(1)
  n_last=np-points_to_discard
  r_last=r(n_last) !last point considering ifin

  ! shift vector rr by one place
  do ir=2,(n_last-1)
     r(ir-1)=r(ir)
  enddo

  ! this is the final grid over rr
  ! r_start--r(1)--r(2)--...--r(n_last-2)--r_last
  !          ^-----------------------^
  !            this the interval passed 
  !              to spline subroutines

  if (debug>=2)then
     write(*,999)'+-------------------------+'
     write(*,999)'|extrapolation and fitting|'
     write(*,999)'+-------------------------+'
  endif
  
  lambda_cicle: do l=1,n_lambda
     outer_vib_level: do n1=1,num_vib_levels
        inner_vib_level:do n2=n1,num_vib_levels
           if (debug>=2)then
              write(*,*)' '
              write(*,1021)'*-->  For symm= ',l,' n= ',n1-1,' np= ',n2-1,'<--*'
           endif
           !------------------------------------------------
           ! *for short range*
           !  v=a/r + b*r
           v_start=v_lambda(1, l, n1, n2)
           v_last =v_lambda(n_last, l, n1, n2)
           ! create v_temp vector for spline, shifting vlam by one place      
           do ir=2,(n_last-1)
              v_temp(ir-1)=v_lambda(ir,l,n1,n2)
              v_lambda(ir-1,l,n1,n2)=v_temp(ir-1)
           enddo
           short_range_derivative=(v_temp(1) - v_start) / (r(1) - r_start)
           a=0.5_wp * r(1) * (v_temp(1) - r(1) * short_range_derivative)
           b=0.5_wp / r(1) * (v_temp(1) + r(1) * short_range_derivative)
           if (debug>=2)then
              write(*,999)'short range parameters v= a/r + b*r '
              write(*,1017)'a = ',  a
              write(*,1018)'b = ',  b
              write(*,*)'---------------------------------------------- '
           endif
                      
           !-----------------------------------------------
           ! *now long range*
           !  -> v=c*r**(-n) + d*r**(-m)
           set_to_zero2: if(l>lambda_to_zero)then
              c=0.0
              d=0.0
              v_lambda(n_last-2, l, n1, n2)=0.0
              v_temp(n_last-2)=0.0
              long_range_derivative=0.0
              cycle  inner_vib_level
           endif set_to_zero2
           long_range_derivative=(v_last - v_temp(n_last-2)) & 
                / (r_last - r(n_last-2))
           
           c=(r(n_last-2) ** n(l)) *  (v_temp(n_last-2) + &
                r(n_last-2) * long_range_derivative / real(m(l), kind=wp))
           
           d=(r(n_last-2) ** m(l)) * (v_temp(n_last-2) + &
                r(n_last-2) * long_range_derivative / real( n(l) ,kind=wp))
           
           d=real(n(l),kind=wp) / real(n(l) - m(l), kind=wp) * d
           c=real(m(l),kind=wp) / real(m(l) - n(l), kind=wp) * c
           if (debug>=2)then
              write(*,999)'long range parameters v=c*r**(-n) + d*r**(-m)'
              write(*,1017)'c = ', c ,'n = ',  n(l) 
              write(*,1018)'d = ', d ,'m = ',  m(l) 
           endif
           
           ! some possible error sources for long range fitting follow:
           
           ! for example if c*d < 0 means that the fitting function has
           ! maximum or a minimum for R>0. If d>0 we have a minimum, if
           ! d<0 we have a maximum.
           ! We look for the point in which we have V'=0. And print a 
           ! warning if this point is larger than r_last
           if( (c*d) <= 0.0_wp ) then
              r_modal= ((-m(l)*d) / (n(l)*c)) ** (n(l)-m(l))
              if (r_modal> r_last) then 
                 !we ask for a check only in the interpolation region
                 write(*,999)'Warning: check long range behavior'
                 write(*,999)'non-monotonic long range function'
                 if (d < 0.0_wp) then
                    write(*,*)'there is a maximum in the interpolating'
                    write(*,*)'function for r= ', r_modal
                 endif
                 if (d > 0.0_wp) then
                    write(*,*)'there is a minimum in the interpolating'
                    write(*,*)'function for r= ', r_modal
                 endif
              endif
           endif
           if (debug>=2)write(*,*)'-----------------------------------------------'
           if (debug>=2)write(*,*)'-----------------------------------------------'
           extrapolation_coeff(1,l,n1,n2)=a
           extrapolation_coeff(2,l,n1,n2)=b
           extrapolation_coeff(3,l,n1,n2)=c
           extrapolation_coeff(4,l,n1,n2)=d
           extrapolation_coeff(1,l,n2,n1)=a
           extrapolation_coeff(2,l,n2,n1)=b
           extrapolation_coeff(3,l,n2,n1)=c
           extrapolation_coeff(4,l,n2,n1)=d
           
           ! generate spline coefficients here
           call spline(r,v_temp,n_last-2, &
                short_range_derivative,long_range_derivative,v2_temp)
           spline_coeff(1:n_last-2,l,n1,n2)=v2_temp(1:n_last-2)
        enddo inner_vib_level
     enddo outer_vib_level
  enddo lambda_cicle
1016 format(t3,a,i3)
1018 format(t3,a,f20.10,5x,a,i3)
1017 format(t3,a,e20.10e3,5x,a,i3)      
1021 format(t3,a,i3,5x,a,i3,5x,a,i3,2x,a)
999 format(t3,a) 
  return
end subroutine vfitt










