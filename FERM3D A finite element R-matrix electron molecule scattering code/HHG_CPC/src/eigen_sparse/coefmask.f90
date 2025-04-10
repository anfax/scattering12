!!!*************************************************************
! 文件/File: coefmask.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: coefmask.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:43
!*************************************************************


SUBROUTINE coefmask(cmask,x,nlftmsk,nrgtmsk)
  USE FiniteElements
  implicit none

  INTEGER :: nlftmsk, nrgtmsk
  TYPE(FEBnodes), INTENT(in) :: x
  REAL(kind=8) :: cmask(*)
  INTEGER :: nbasis, nmnds, j
  REAL(kind=8), ALLOCATABLE :: ndmsk(:)

  nmnds = x%numnds
  nbasis = x%numbasis
  ALLOCATE( ndmsk(nmnds) )
  ndmsk(:) = 1.d0
  cmask(1:nbasis) = 1.d0

  DO j = 1,nlftmsk
    ndmsk(j) = nodemask(x%nodes(j),x%nodes(nlftmsk),x%nodes(1))
  END DO

  DO j = (nmnds-nrgtmsk+1),nmnds
    ndmsk(j) = nodemask(x%nodes(j),x%nodes(nmnds-nrgtmsk+1),x%nodes(nmnds))
  END DO

  WRITE(*,*)
  WRITE(*,*) 'nodemask:'
  DO j = 1,nlftmsk+1
    WRITE(*,*) x%nodes(j), ndmsk(j)
  END DO
  DO j = nmnds-nrgtmsk,nmnds
    WRITE(*,*) x%nodes(j), ndmsk(j)
  END DO
  WRITE(*,*)

  cmask(1:x%numleft) = ndmsk(1)
  DO j = 2,(nmnds-1)
    cmask((3*j-5+x%numleft):(3*j-3+x%numleft)) = ndmsk(j)
  END DO
  cmask((nbasis-x%numright+1):nbasis) = ndmsk(nmnds)

  RETURN

  CONTAINS

  REAL(kind=8) FUNCTION nodemask(xx,ndi,ndx)
    REAL(kind=8) :: xx, ndi, ndx, zx, prec, ac, dz

    prec = 0.01d0
    dz = abs(ndx-ndi)
    IF ( dz < 1.d-40 ) THEN
      nodemask = 1.d0
    ELSE
      zx = abs(ndx-xx)
      ac = -2.d0*log(2.d0*prec)/log(2.d0-2.d0*prec)/dz 
      nodemask = 0.5d0*(1.d0+ tanh(ac*(zx-0.5d0*dz))) 
    ENDIF

    RETURN
  END FUNCTION nodemask


  REAL(kind=8) FUNCTION old_nodemask(node,ndi,ndx)
    REAL(kind=8) :: node, ndi, ndx, sclf

    sclf = 4.d0
    old_nodemask = tanh(sclf*(node-ndx)/(ndi-ndx))
    old_nodemask = (old_nodemask - tanh(0.d0))/tanh(sclf)
    RETURN
  END FUNCTION old_nodemask

END SUBROUTINE coefmask
