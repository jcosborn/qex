#!/bin/sh

c2n="c2nim"
#c2n="/home/josborn/work/lqcd/src/c2nim-mygit/c2nim"
#QUDA=$HOME/lqcd/src/quda-0.8.0
QUDA=$HOME/work/lqcd/src/quda-git
QI=$QUDA/include

ezero=`grep -A1 enum $QI/enum_quda.h |grep -v enum |grep -v = |grep QUDA |sed 's/,.*//' |sort`
#echo $ezero
szero=`echo "$ezero" |sed 's|\(QUDA_[A-Z0-9_]*\)|s/\1/\1=0/;|'`
echo $szero

cat $QI/enum_quda.h |\
    sed 's/= *QUDA_INVALID_ENUM/= -2147483647/' >enum_quda.cnim
$c2n --header -o:enum_quda.cnim2 quda.c2nim enum_quda.cnim
cat enum_quda.cnim2 |\
    sed 's/= -2147483647/= QUDA_INVALID_ENUM/' |
    sed "$szero" >enum_quda_new.nim
rm enum_quda.cnim enum_quda.cnim2

$c2n --header -o:quda_constants_new.nim quda.c2nim $QI/quda_constants.h

#cat $QI/quda.h |sed 's/double _Complex/dcomplex/' >quda.cnim
cat $QI/quda.h |sed '/#ifndef __CUDACC_RTC__/,/#endif/c\
typedef double double_complex[2];\
' >quda.cnim
$c2n --header -o:quda.cnim2 quda.c2nim quda.cnim
cat <<EOF >quda_new.nim
import enum_quda, quda_constants
EOF
cat quda.cnim2 >>quda_new.nim
rm quda.cnim quda.cnim2


exit

cat enum_quda_fortran.h |sed 's|^!|//|' | \
 sed 's/QudaMassNormalization/QudaMassNormalizationType/' | \
 sed 's/^#define  *\(Quda[^ ]*\)  *integer(4) *$/typedef int \1;/' | \
 sed 's/^#define  *\(Quda[^ ]*\) *$/typedef int \1;/' >enum_quda.cnim
c2nim --header enum_quda.cnim

#cat enum_quda.h |sed 's/INT_MIN/low(cint)/' >enum_quda.cnim
#c2nim --header enum_quda.cnim

c2nim --header quda_constants.h

cat quda.h |sed '
 s/#include <enum_quda.h>/#@import enum_quda\n@#/
 s/#include <quda_constants.h>/#@import quda_constants\n@#/
 s/QudaMassNormalization/QudaMassNormalizationType/
' >quda_main.cnim
c2nim --header quda_main.cnim
