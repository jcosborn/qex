c2n=${HOME}/.nimble/bin/c2nim
sed_inplace=sed -i
# bsd sed needs an (empty string) argument to -i
#sed_inplace=sed -i ''
#QI=../../../quda/quda/include
#QI=../../../quda-git/include
QI=../../../quda-git-gk/include

test:\
	enum_quda.nim.test\
	quda.nim.test\
	quda_constants.nim.test\
	quda_milc_interface.nim.test

all:\
	enum_quda.nim\
	quda.nim\
	quda_constants.nim\
	quda_milc_interface.nim\

# doesn't work
#	color_spinor_field.nim\
#	host_utils.nim\

.PHONY:	all

%.nim: %.nim.test
	mv -i $@ $@.orig
	mv $< $@

enum_quda.nim.test: quda.c2nim ${QI}/enum_quda.h
	cat $^ |sed -e 's/= *QUDA_INVALID_ENUM/= -2147483647/'\
		-e '/[^a-zA-Z0-9_,;=]enum[^a-zA-Z0-9_,;=]/{:j;/[,;]/!{N;bj;};s/\([^a-zA-Z0-9_,;=]enum[^a-zA-Z0-9_,;=][^,;=]*\),/\1 = 0,/;}'\
		>enum_quda.h.tmp
	${c2n} -o:$@ --header enum_quda.h.tmp
	${sed_inplace} -e 's/enum_quda\.h\.tmp/enum_quda.h/g'\
		-e 's/-2147483647/QUDA_INVALID_ENUM/g'\
		-e '/^const$$/{:j;/\n$$/!{N;bj;};s/\([a-zA-Z0-9]\) =/\1\* =/g;}'\
		$@
	rm enum_quda.h.tmp

quda_constants.nim.test: quda.c2nim ${QI}/quda_constants.h
	cat $^ >quda_constants.h
	${c2n} -o:$@ --header quda_constants.h
	rm quda_constants.h

quda_milc_interface.nim.test: quda.c2nim ${QI}/quda_milc_interface.h
	cat $^ |sed '/#if __COMPUTE_CAPABILITY__/,/#endif/d' >quda_milc_interface.h
	${c2n} -o:$@ --header quda_milc_interface.h
	${sed_inplace} -e '1i\
	import enum_quda, quda\
	' $@
	echo >>$@
	rm quda_milc_interface.h

quda.nim.test: quda.c2nim ${QI}/quda.h
	cat $^ |sed -e 's/const int/ConstInt/g' \
	-e '/#ifndef __CUDACC_RTC__/,/#endif/c\
	#@\
	type\
	  ConstInt* {.importc:"const int".} = cint\
	  double_complex* {.importc:"double _Complex".} = object\
	converter toDoubleComplex*(x: array[2,float]): double_complex =\
	  var r = cast[ptr array[2,float]](addr result)\
	  r[] = x\
	@#\
	' >quda.h.tmp
	${c2n} -o:$@ --header quda.h.tmp
	${sed_inplace} -e '1i\
	import enum_quda, quda_constants\
	'\
		-e '/^type$$/{:j;/\n$$/!{N;bj;};s/\(= *proc[^a-zA-Z0-9_\n]*[^}]*\)\n/\1 {.cdecl.}\n/g;}'\
		-e 's/quda\.h\.tmp/quda.h/g'\
		$@
	rm quda.h.tmp

host_utils.nim:	quda.c2nim ${QI}/../tests/utils/host_utils.h
	sed 's,^\(extern [a-zA-Z0-9_]*\) & *\([a-zA-Z0-9_]*;\),\1 \2,' \
		${QI}/../tests/utils/host_utils.h > host_utils.h.tmp
	${c2n} -o:$@ --header quda.c2nim host_utils.h.tmp
	${sed_inplace} -e '1i\
	import quda, color_spinor_field\
	'\
		-e 's/host_utils\.h\.tmp/host_utils.h/g'\
		host_utils.nim
	rm host_utils.h.tmp

color_spinor_field.nim:	quda.c2nim ${QI}/color_spinor_field.h
	${c2n} -o:$@ --header $^
