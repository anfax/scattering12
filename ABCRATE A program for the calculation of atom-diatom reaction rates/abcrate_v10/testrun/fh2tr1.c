/**************************************************************
 * 文件/File: fh2tr1.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-11 00:14:06
 **************************************************************
 */

/**************************************************************
 * 文件/File: fh2tr1.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-10 23:45:42
 **************************************************************
 */

#!/bin/csh -f
#
#usage fh2tr1.c
#
cd ..
set root   = `pwd`
#
set abcdir = $root/src
set exedir = $root/exe
set potdir = $root/poten
set wrkdir = $root/testrun
set scrdir = $root/script
#
set surf    = fh2m5
set abcrate = $surf.abc
set ofile   = fh2tr1
#
set CF = `$scrdir/comp.c compile`
set LD = `$scrdir/comp.c load`
if ('$CF' == 'error' || '$LD' == 'error') then
  echo "Unable to determine proper compiler/loader commands."
  exit
endif
#
cd $wrkdir
#
# Check to see if the executable file exist, if not make it.
#
if (! -f $exedir/$surf.abc) then
    if (! -f $potdir/$surf.o) then
        cd $potdir
        $CF $surf.f
    endif
    if (! -f $potdir/abc.o) then
        cd $potdir
        $CF abc.f
    endif
    cd $wrkdir
    $LD $exedir/$surf.abc $potdir/$surf.o $potdir/abc.o $abcdir/*.o 
endif
#
cp fh2tr1.dat abc.5
#
time $exedir/$abcrate 
#
# rename abcrate output files
#
if (-e abc.5)  rm abc.5
if (-e abc.6)  mv abc.6  $ofile.out
if (-e abc.1)  mv abc.1  $ofile.1
if (-e abc.11) mv abc.11 $ofile.11
if (-e abc.12) mv abc.12 $ofile.12
if (-e abc.14) mv abc.14 $ofile.14
if (-e abc.20) mv abc.20 $ofile.20
if (-e abc.21) mv abc.21 $ofile.21
if (-e abc.22) mv abc.22 $ofile.22
#
