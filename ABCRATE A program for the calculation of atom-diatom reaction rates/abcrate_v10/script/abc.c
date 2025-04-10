/**************************************************************
 * 文件/File: abc.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-11 00:14:06
 **************************************************************
 */

/**************************************************************
 * 文件/File: abc.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-10 23:45:40
 **************************************************************
 */

#!/bin/csh -f
#
# If the user invokes this script with "abc.c clean" then use the "clean"
# target of the abcmake file.
if ($1 == "clean") then
  cd ../src
  make -f abcmake clean
  cd ../poten
  rm -f *.o
  exit
endif
#
# Get compiler command and machine type, set environment variables for the
# makefile, change directory to the source directory, and execute the makefile.
#
set cf   = `comp.c compile`
set mach = `comp.c machine`

setenv CF   "$cf"
setenv MACH "$mach"

cd ../src

make -f abcmake
#
