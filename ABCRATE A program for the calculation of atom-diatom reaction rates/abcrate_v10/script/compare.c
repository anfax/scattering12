/**************************************************************
 * 文件/File: compare.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-11 00:14:06
 **************************************************************
 */

/**************************************************************
 * 文件/File: compare.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-10 23:45:40
 **************************************************************
 */

#!/bin/csh -f
#
# Compile compare.f if that hasn't been done yet, make the cmp directory in
# the testrun directory if possible, compare all <testrun>.22 output in
# the testrun directory to that stored in the testo directory.  The output
# of the comparison program is stored in the testrun/cmp directory.
#
set script = `pwd`
#
if (! -e compare.x) then
  set cf = `comp.c compile`
  set ld = `comp.c load`
  $cf compare.f
  $ld compare.x compare.o
  if (-e compare.o) /bin/rm -f compare.o
endif

if (! -e compare.x) then
  echo "error:  can't find compare.x"
  exit
endif
#
cd ..
set root=`pwd`

set dir1 = $root/testo
set dir2 = $root/testrun

set odir = $dir2/cmp
if (-e $odir) then
  echo "error:  comparison output directory already exists."
  exit
endif
mkdir $odir

cd $script

foreach x ($dir2/*.22)
  set fname = $x:t
  if (-e $dir1/$fname) then
    cp $dir1/$fname ./file1
    cp $dir2/$fname ./file2
    set oname = $odir/$fname:r.cmp
    ./compare.x > $oname
    /bin/rm -f file1
    /bin/rm -f file2
  endif
end
#
cd $script
/bin/rm -f compare.x
#
