#!/bin/sh

c2n="c2nim"
#c2n="/home/josborn/work/lqcd/src/c2nim-mygit/c2nim"
QUDA=$HOME/lqcd/src/quda-0.8.0
#QUDA=$HOME/work/lqcd/src/quda-git
QI=$QUDA/include

cat $QI/enum_quda.h |\
    sed 's/= *QUDA_INVALID_ENUM/= -2147483647/' >enum_quda.cnim
$c2n --header -o:enum_quda.cnim2 quda.c2nim enum_quda.cnim
cat enum_quda.cnim2 |\
    sed 's/= -2147483647/= QUDA_INVALID_ENUM/' >enum_quda.nim
rm enum_quda.cnim enum_quda.cnim2

$c2n --header -o:quda_constants.nim quda.c2nim $QI/quda_constants.h

$c2n --header -o:quda.cnim quda.c2nim $QI/quda.h
cat <<EOF >quda.nim
import enum_quda, quda_constants
EOF
cat quda.cnim >>quda.nim
rm quda.cnim


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
