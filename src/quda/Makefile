c2n=${HOME}/.nimble/bin/c2nim
sed_inplace=sed -i
# bsd sed needs an (empty string) argument to -i
#sed_inplace=sed -i ''
QI=../../../quda/quda/include

all:\
	enum_quda.nim\
	quda.nim\
	quda_constants.nim\
	quda_milc_interface.nim\

# doesn't work
#	color_spinor_field.nim\
#	host_utils.nim\

.PHONY:	all

enum_quda.nim:	quda.c2nim ${QI}/enum_quda.h
	sed -e 's/= *QUDA_INVALID_ENUM/= -2147483647/'\
		-e '/[^a-zA-Z0-9_,;=]enum[^a-zA-Z0-9_,;=]/{:j;/[,;]/!{N;bj;};s/\([^a-zA-Z0-9_,;=]enum[^a-zA-Z0-9_,;=][^,;=]*\),/\1 = 0,/;}'\
		${QI}/enum_quda.h >enum_quda.h.tmp
	${c2n} -o:$@ --header quda.c2nim enum_quda.h.tmp
	${sed_inplace} -e 's/enum_quda\.h\.tmp/enum_quda.h/g'\
		-e 's/-2147483647/QUDA_INVALID_ENUM/g'\
		-e '/^const$$/{:j;/\n$$/!{N;bj;};s/\([a-zA-Z0-9]\) =/\1\* =/g;}'\
		enum_quda.nim
	rm enum_quda.h.tmp

quda_constants.nim:	quda.c2nim ${QI}/quda_constants.h
	${c2n} -o:$@ --header $^

quda_milc_interface.nim: quda.c2nim ${QI}/quda_milc_interface.h
	${c2n} -o:$@ --header $^
	${sed_inplace} -e '1i\
	import enum_quda, quda\
	' quda_milc_interface.nim

quda.nim:	quda.c2nim ${QI}/quda.h
	sed '/#ifndef __CUDACC_RTC__/,/#endif/c\
	typedef double double_complex[2];\
	' ${QI}/quda.h >quda.h.tmp
	${c2n} -o:$@ --header quda.c2nim quda.h.tmp
	${sed_inplace} -e '1i\
	import enum_quda, quda_constants\
	'\
		-e '/^type$$/{:j;/\n$$/!{N;bj;};s/\(= *proc[^a-zA-Z0-9_\n]*[^}]*\)\n/\1 {.cdecl.}\n/g;}'\
		-e 's/quda\.h\.tmp/quda.h/g'\
		quda.nim
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
