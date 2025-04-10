!!!*************************************************************
! 文件/File: eigenstuff.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: eigenstuff.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


SUBROUTINE ConvergenceInfo(info, nbasis, nev, ncv, which, nconv, iparam, tol)
  implicit none

  INTEGER :: info, nbasis, nev, ncv, nconv
  INTEGER :: iparam(11)
  CHARACTER :: which*2
  REAL(kind=8) :: tol

!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%

  IF ( info == 1) THEN
    print *, ' '
    print *, ' Maximum number of iterations reached.'
    print *, ' '
  ELSE IF ( info == 3) THEN
    print *, ' ' 
    print *, ' No shifts could be applied during implicit', &
          & ' Arnoldi update, try increasing NCV.'
    print *, ' '
  END IF      

  print *, ' '
  print *, ' Convergence Information'
  print *, ' ======================= '
  print *, ' '
  print *, ' Size of the matrix is ', nbasis
  print *, ' The number of Ritz values requested is ', nev
  print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
  print *, ' What portion of the spectrum: ', which
  print *, ' The number of converged Ritz values is ', nconv 
  print *, ' The number of Implicit Arnoldi update', &
        & ' iterations taken is ', iparam(3)
  print *, ' The number of OP*x is ', iparam(9)
  print *, ' The convergence criterion is ', tol
  print *, ' '

 RETURN
END SUBROUTINE ConvergenceInfo



SUBROUTINE sorteigs(idx,eig,nev)
  implicit none

  INTEGER :: nev, idx(nev)
  REAL(kind=8) :: eig(nev)
  INTEGER :: j, itmp, count
  REAL(kind=8) :: sofar

  DO j = 1,nev
    idx(j) = j
  END DO

  DO count = 2,nev
    DO j = 2,nev
      IF (eig(idx(j)) < eig(idx(j-1))) THEN
        itmp = idx(j)
        idx(j) = idx(j-1)
        idx(j-1) = itmp
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE sorteigs
