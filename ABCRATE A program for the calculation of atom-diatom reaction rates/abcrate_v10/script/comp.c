/**************************************************************
 * 文件/File: comp.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-11 00:14:06
 **************************************************************
 */

/**************************************************************
 * 文件/File: comp.c
 * 由李皓(lihao,LIH_ao@outlook.com)创建/修改/分享
 * Created/Modified/Shared by Li Hao (lihao,LIH_ao@outlook.com)
 * 日期时间/Date Time: 2025-04-10 23:45:40
 **************************************************************
 */

#!/bin/csh -f
#
# args:
#   $1 - either "compile" to return the appropriate compiler commands
#            or "load"    to return the appropraite loader command
#            or "machine" to return the machine type
#   $2 - architecture/machine type (optional)
#        this may be useful for unsupported architectures or custom compiles
#
# This is a simple script which determines the type of architecture it is
# running on.  It uses the UNIX commands "uname" and "nawk" or "awk".  The
# script is designed to detect Cray, IBM RS/6000, and SGI Iris machines.
# This script will need to be modified for other architectures.
#
# Let's use "nawk" if it is present on this system, otherwise try "awk" if
# available.  If neither is available, exit with an error.
# We will look for nakw/awk in /bin and /usr/bin.  If it's not there, we quit,
# unless an architecture has been specified.
#
if ($2 != "") then
  set arch = $2
  set mach = $2
else
  if (-e /bin/nawk) then
    set awk = /bin/nawk
  else if (-e /usr/bin/nawk) then
    set awk = /usr/bin/nawk
  else if (-e /bin/awk) then
    set awk = /bin/awk
  else if (-e /usr/bin/awk) then
    set awk = /usr/bin/awk
  else
    echo "error"
    exit
  endif
#
# We will get the architecture information from the "uname" command.  We use
# nawk/awk to parse the string we get back and find the architecture info.
#
  set arch = ` uname -a | $awk ' { if (/IRIX/ || /AIX/) {print $1} else { if (/CRAY/) {print $5} } } ' `
#
# Now process the architecture information and return "sgi", "ibm", or "cray".
# If none of these machines are successfully detected and the user did not
# enter another name, exit with an error.
#
  if ($arch == 'IRIX' || $arch == 'IRIX64') then
    set mach = 'sgi'
  else if ($arch == 'AIX') then
    set mach = 'ibm'
  else if ($arch == 'CRAY') then
    set mach = 'cray'
  else if ($mach == '')
    echo "error"
    exit
  endif
endif
#
# Finally, give the compiler and loader commands appropriate for the selected
# architecture.  If the machine type is not recognized, exit with an error.
#
if ($mach == 'sgi') then
  set cmp = "f77 -c -64 -mips4 -O3 -static"
  set ldr = "f77 -64 -o "
else if ($mach == 'ibm') then
  set cmp = "xlf -c -O3 -qstrict -qhot -qsave"
  set ldr = "xlf -o "
else if ($mach == 'cray') then
  set cmp = "cft77 -dp -astatic -i64"
  set ldr = "segldr -o "
else
  echo "error"
  exit
endif
#
if ($1 == 'compile') then
  echo $cmp
else if ($1 == 'load') then
  echo $ldr
else if ($1 == 'machine') then
  echo $mach
else
  echo "error"
endif
#
