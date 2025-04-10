!!!*************************************************************
! 文件/File: makebessel.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: makebessel.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************



SUBROUTINE makebessel(bessj, x, mxn)
  implicit none

  INTEGER :: mxn
  REAL(kind=8) :: bessj(0:mxn), x
  INTEGER :: n, j, mxtrms, ntrms, mu
  REAL(kind=8), ALLOCATABLE :: terms(:)
  REAL(kind=8) :: tmp, pi, tol

  pi = 4.d0*atan(1.d0)
  tol = 1.d-30
  mxtrms = mxn + INT(SQRT(.5d0*mxn))
  mxtrms = MAX(100, mxtrms)
  ALLOCATE( terms(0:mxtrms) )
  terms(:) = 0.d0

  IF ( abs(x) < 0.35d0 ) THEN
    !WRITE(*,*) 'using polynomial approximation for bessel functions'
    DO n = 0,mxn
      IF (n==0) THEN
        tmp = 1.d0
      ELSE
        tmp = 0.5d0*x*tmp/n
      ENDIF
      bessj(n) = 1.d0 - x*x / (12.d0*(n+3))
      bessj(n) = 1.d0 - x*x*bessj(n) / (8.d0*(n+2))
      bessj(n) = 1.d0 - x*x*bessj(n) / (4.d0*(n+1))
      bessj(n) = tmp*bessj(n)
    END DO
  ELSE
    !WRITE(*,*) 'using recurrence relation for bessel functions'
    IF ( mxn < 500 ) THEN
      tmp = 0.5d0*abs(x)
      terms(0) = 1.d0
      terms(1) = tmp
      DO j = 2,mxn
        terms(j) = tmp*terms(j-1) / j
      END DO
      ntrms = mxn+1
      tmp = 0.5d0*x
      DO j = ntrms,mxtrms
        terms(j) = tmp*terms(j-1) / j
        IF ( terms(j) < tol ) THEN
          ntrms = j
          EXIT
        ENDIF
      END DO
    ELSE
      ntrms = mxtrms
      terms(ntrms) = 1.d-100
    ENDIF

    if ( (x<0.d0).AND.(MOD(ntrms,2)==1) ) terms(ntrms) = -terms(ntrms)

    terms(ntrms-1) = terms(ntrms)*(2*ntrms/x)
    DO n = (ntrms-1),1,-1
      terms(n-1) = terms(n)*(2*n/x) - terms(n+1)
    END DO
    tmp = 0.d0
    DO j = (ntrms-MOD(ntrms,2)),2,-2
      tmp = tmp + 2.d0*terms(j)
    END DO
    tmp = tmp + terms(0)
    tmp = 1.d0/tmp
    terms(0:ntrms) = tmp*terms(0:ntrms)
    terms(ntrms+1:mxtrms) = 0.d0
    bessj(0:mxn) = terms(0:mxn)
  ENDIF

  DEALLOCATE( terms )
  RETURN
END SUBROUTINE makebessel
