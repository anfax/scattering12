!!!*************************************************************
! 文件/File: smallF.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: smallF.f90
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:45
!*************************************************************

module smallFunctions
contains

! Function for factorial, Used by PLGNDR Function   
    integer*8 function fact(x)
    integer*8 x,i

    if (x.eq.0) then
        fact=1
        return
    else if (x.eq.1) then
        fact=1
        return
    else
        fact=1
        do i=1,x
        fact=fact*i
        end do
    end if
    
    return
    end function fact

! Function for Double factorial 
    integer*8 function dfact(x)
    integer*8 x,i

    if (x.eq.1) then
        dfact=1
        return
    else
        dfact=1
        do i=1,x,2
        dfact=dfact*i
        end do
    end if

    return
    end function

end module
