!!!*************************************************************
! 文件/File: dipele_general.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-11 00:14:07
!*************************************************************

!!!*************************************************************
! 文件/File: dipele_general.f
! 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
! Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
! 日期时间/Date Time: 2025-04-10 23:45:44
!*************************************************************

       subroutine dipele(r1,r2,cgam,dx,dy,dz,nelec,nelecmax)
       implicit real*8(a-h,o-z)

       ! r1 is the distance 01 (small r in Jacobi coordinates)
       ! r2 is the distance between the c-o-m of 01 to atom 2 (big R in
       !                                        Jacobi coordinates
       ! cgam is cos(gamma), gamma being the angle between r1 and r2
       !                      vectors
       ! distances and dipole are in atomic units

       dimension dx(nelecmax),dy(nelecmax),dz(nelecmax)

        return
        end

        subroutine setdipini
        implicit real*8(a-h,o-z)

        return
        end
