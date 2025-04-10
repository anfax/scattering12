!!!*************************************************************
! 文件/File: vpot.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:06
!*************************************************************

!!!*************************************************************
! 文件/File: vpot.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:42
!*************************************************************

SUBROUTINE vpot(rp, n, u, diag)
  ! this subroutine returns u: the matrix potential
  ! and diag: the diagonal part of it in different vectors
  USE PRECISION
  USE shared
  ! rename num_pot_point to npun and num_lambda_terms to n_lambda
  ! rename extrapolation_coeff to ce
  USE potential, npun=>num_pot_points, n_lambda=>num_lambda_terms, &
       ce => extrapolation_coeff
  use long_range_mod
  use percival_seaton, coeff=>pot_coeff  !potential coefficents
  use channels
  use potential_scratch ! this 'use' give access to a scratch area 
                        ! that contains vtt v_tmp and v_tmp2
                        ! these vectors are of constant size and
                        ! thus allocated in main() 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(kind=wp), INTENT(IN) :: rp
  REAL(kind=wp), INTENT(OUT), DIMENSION(n,n) :: u
  REAL(kind=wp), INTENT(OUT), DIMENSION(n) :: diag

  REAL(kind=wp) :: over_r, vl, sum, x
  INTEGER :: n_last, exp, exp_2, n1, n2, l, i, ii, nmax, myrank
  !    +-----------------------------------------------+
  !    | note (some differences with molscat):         |  
  !    | k^2 matrix is handled inside propagator.      |
  !    | molscat sum also the difference eint-ered but |
  !    | both quantities are setted to 0.0 by daprop.  |
  !    | in our code we prefer put k^2 diagonal matrix |
  !    | directly inside daprop.                       |
  !    +-----------------------------------------------+
  !     this subroutine must be fast so we play with indexes
  !     and we avoid every write statement.
  nmax=vib_levels
  over_r=1.0_wp/rp/rp
  n_last= npun - points_to_discard 
  !     vlambdas are mapped only by vibrational levels.
  !     so calculate {v(r)} only for the couples n n'
  !     this reduce overhead in pow_and_exp built in 
  !     fortran procedure.
  
  IF ( rp > r(n_last-2) ) THEN   ! long range
     DO n1=1,nmax
        DO n2=n1,nmax
           DO l=1,n_lambda
              vtt(l,n1,n2)=ce(3,l,n1,n2)/(rp**n_exp(l)) + &
                   ce(4,l,n1,n2)/(rp**m_exp(l))
              vtt(l,n2,n1)=vtt(l,n1,n2)
           ENDDO
        ENDDO
     ENDDO
  ELSEIF(rp<r(1)) THEN    ! short range
     DO n1=1,nmax
        DO n2=n1,nmax
           DO l=1,n_lambda
              vtt(l,n1,n2)=ce(1,l,n1,n2)/rp + ce(2,l,n1,n2)*rp
              vtt(l,n2,n1)=vtt(l,n1,n2)
           ENDDO
        ENDDO
     ENDDO
  ELSE                   ! intermediate range
     DO n1=1,nmax
        DO n2=n1,nmax
           DO l=1,n_lambda
              v_tmp (1:n_last-2)=v_lambda(1:n_last-2,l,n1,n2)
              v_tmp2(1:n_last-2)=spline_coeff(1:n_last-2,l,n1,n2)
              CALL splint(r, v_tmp, v_tmp2, n_last-2, rp, vl)
              vtt(l,n1,n2)=vl
              vtt(l,n2,n1)=vl
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  !     calculate couplings
  DO ii=1,n
     DO i=ii,n
        sum=0.0d0
        DO l=1,n_lambda
           sum=sum+coeff(l,i,ii)*vtt( l, n_local(i)+1, n_local(ii)+1 )
        ENDDO
        x=two_mu*sum
        IF (i == ii)THEN
           x = x + l_local(i)*(l_local(i)+1) * over_r
           diag(i)=x
        ENDIF
        u(i,ii)=x
        u(ii,i)=x
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE vpot






