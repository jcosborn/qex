import simdWrap
export simdWrap

{.pragma: imm, header:"immintrin.h".}
{.pragma: imms, header:"immintrin.h", incompleteStruct.}
type
  m64*   {.importc: "__m64"  , imms.} = object
  m128*  {.importc: "__m128" , imms.} = object
  m128d* {.importc: "__m128d", imms.} = object
  m128i* {.importc: "__m128i", imms.} = object
  m128h* = distinct int64
  m256*  {.importc: "__m256" , imms.} = object
  m256d* {.importc: "__m256d", imms.} = object
  m256i* {.importc: "__m256i", imms.} = object
  m256h* = distinct m128i
  m512*  {.importc: "__m512" , imms.} = object
  m512d* {.importc: "__m512d", imms.} = object
  m512i* {.importc: "__m512i", imms.} = object
  m512h* = distinct m256i
  mmask8*  {.importc: "__mmask8" , imms.} = object
  mmask16* {.importc: "__mmask16", imms.} = object
  mmask32* {.importc: "__mmask32", imms.} = object
  mmask64* {.importc: "__mmask64", imms.} = object
#[
{.pragma: imm, header:"immintrin.h", incompleteStruct.}
type
  m64*   {.importc: "__m64"  , imm.} = object
    a: array[2,float32]
  m128*  {.importc: "__m128" , imm.} = distinct array[4,float32]
  m128d* {.importc: "__m128d", imm.} = distinct array[2,float64]
  m128i* {.importc: "__m128i", imm.} = object
    a: array[4,int32]
  m128h* = distinct int64
  m256*  {.importc: "__m256" , imm.} = object
    a: array[8,float32]
  m256d* {.importc: "__m256d", imm.} = distinct array[4,float64]
  m256i* {.importc: "__m256i", imm.} = object
    a: array[8,int32]
  m256h* = distinct m128i
  m512*  {.importc: "__m512" , imm.} = distinct array[16,float32]
  m512d* {.importc: "__m512d", imm.} = distinct array[8,float64]
  m512i* {.importc: "__m512i", imm.} = distinct array[16,int32]
  m512h* = distinct m256i
  mmask8*  {.importc: "__mmask8" , imm.} = object
  mmask16* {.importc: "__mmask16", imm.} = object
  mmask32* {.importc: "__mmask32", imm.} = object
  mmask64* {.importc: "__mmask64", imm.} = object
]#

type
  SimdX86S* = m64 | m128 | m256 | m512
  SimdX86D* = m128d | m256d | m512d
  SimdX86* = SimdX86S | SimdX86D

when defined(SSE):
  type
    SimdS4* = Simd[m128]
    SimdD2* = Simd[m128d]
    SimdI4* = Simd[m128i]
    SimdH4* = Simd[m128h]
when defined(AVX):
  type
    SimdS8* = Simd[m256]
    SimdD4* = Simd[m256d]
    SimdI8* = Simd[m256i]
    SimdH8* = Simd[m256h]
when defined(AVX512):
  type
    SimdS16* = Simd[m512]
    SimdD8*  = Simd[m512d]
    SimdI16* = Simd[m512i]
    SimdH16* = Simd[m512h]

template eval*(x: SimdX86): untyped = x

#var CMP_EQ_OS {.importc: "_CMP_EQ_OS", imm.} = cint
var CMP_LT_OS* {.importc: "_CMP_LT_OS", imm.}: cint
