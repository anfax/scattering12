/**************************************************************
 * 文件/File: runtest.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-11 00:14:06
 **************************************************************
 */

/**************************************************************
 * 文件/File: runtest.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-10 23:45:40
 **************************************************************
 */

#!/bin/csh -f
#
#  Execute all test runs is $1 == 'all' or a specific test run if $1 == [1-8]
#  Print out a list of options if $1 is blank.
#
if ($1 == 'all') then
  echo "running all test runs"
  cd ../testrun
  foreach x (*.c)
    echo "running test run " $x:r
    ./$x
  end
  echo "test runs completed"
  exit
else if ($1 == '1') then
  set tname = clhbrtr1.c
else if ($1 == '2') then
  set tname = clhcltr1.c
else if ($1 == '3') then
  set tname = clhcltr2.c
else if ($1 == '4') then
  set tname = dfhtr1.c
else if ($1 == '5') then
  set tname = hfhtr1.c
else if ($1 == '6') then
  set tname = dh2tr1.c
else if ($1 == '7') then
  set tname = dh2tr2.c
else if ($1 == '8') then
  set tname = fh2tr1.c
else if ($1 == '9') then
  set tname = ih2tr1.c
else if ($1 == '10') then
  set tname = oh2tr1.c
else
  echo "please input a test number (1-9) or all"
  echo "  1)  clhbr tr1"
  echo "  2)  clhcl tr1"
  echo "  3)  clhcl tr2"
  echo "  4)  dfh   tr1"
  echo "  5)  hfh   tr1"
  echo "  6)  dh2   tr1"
  echo "  7)  dh2   tr2"
  echo "  8)  fh2   tr1"
  echo "  9)  ih2   tr1"
  echo " 10)  oh2   tr1"
  echo "  all - runs all test runs"
  exit
endif

if (-e ../testrun/$tname) then
  echo "now running test run " $tname:r "..."
  cd ../testrun
  ./$tname
  echo "...test run completed"
else
  echo "error:  unable to find testrun"
endif
#
