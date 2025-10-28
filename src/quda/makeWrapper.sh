#!/bin/sh

#c2n="c2nim"
c2n=/g/g91/jin8/.nimble/bin/c2nim
#c2n="/home/josborn/work/lqcd/src/c2nim-mygit/c2nim"
#QUDA=$HOME/lqcd/src/quda-0.8.0
#QUDA=$HOME/work/lqcd/src/quda-git
QUDA=/p/lustre5/jin8/soft.20251023.rocm6.4.2/quda
QI=$QUDA/include

ezero=`grep -A1 enum $QI/enum_quda.h |grep -v '^ *typedef.*enum' |sed 's,\(//\|/\*\).*,,' |grep -v = |sed -n '/QUDA/{s/^ *//;s/,.*//;p}' |sort`
ezero1=`sed -n '/enum.*{.*}.*;/{s/.*typedef.*enum.*{ *//; s/ *[,}].*//; s,\(//\|/\*\).*,,; p}' $QI/enum_quda.h |grep -v = |sort`
#echo "$ezero"
#echo "$ezero1"
szero=`{ echo "$ezero" ; echo "$ezero1" ; } |sed 's|\(QUDA_[A-Z0-9_]*\)|s/\1\\\([^_a-zA-Z0-9]\\\)/\1 = 0\\\\1/;|'`
echo $szero

cat $QI/enum_quda.h |\
    sed 's/= *QUDA_INVALID_ENUM/= -2147483647/' >enum_quda.cnim
$c2n --header -o:enum_quda.cnim2 quda.c2nim enum_quda.cnim
cat enum_quda.cnim2 |\
    sed 's/= -2147483647/= QUDA_INVALID_ENUM/' |
    sed "$szero" >enum_quda_new.nim
rm enum_quda.cnim enum_quda.cnim2
echo 'WARNING:
	check QudaUpdateSplitGauge in enum_quda.nim
	make sure to have
	const QUDA_UPDATE_SPLIT_GAUGE_OFF = QUDA_UPDATE_SPLIT_GAUGE_FALSE
	or something the same as in enum_quda.h'

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
echo 'WARNING:
	check quda.nim
	make sure to have
	type
	  ConstInt* {.importc:"const int".} = cint
	  double_complex* {.importc:"double _Complex".} = object
	converter toDoubleComplex*(x: array[2,float]): double_complex =
	  var r = cast[ptr array[2,float]](addr result)
	  r[] = x
	type
	  QudaCommsMap* = proc (coords: ptr cint; fdata: pointer): cint {.cdecl.}'


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
