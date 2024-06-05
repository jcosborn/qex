import simdWrap
export simdWrap

{.pragma: imm, header:"immintrin.h".}
type
  m64*   {.importc: "__m64"  , imm.} = object
  m128*  {.importc: "__m128" , imm.} = object
  m128d* {.importc: "__m128d", imm.} = object
  m128i* {.importc: "__m128i", imm.} = object
  m128h* = distinct int64
  m256*  {.importc: "__m256" , imm.} = object
  m256d* {.importc: "__m256d", imm.} = object
  m256i* {.importc: "__m256i", imm.} = object
  m256h* = distinct m128i
  m512*  {.importc: "__m512" , imm.} = object
  m512d* {.importc: "__m512d", imm.} = object
  m512i* {.importc: "__m512i", imm.} = object
  m512h* = distinct m256i
  mmask8*  {.importc: "__mmask8" , imm.} = object
  mmask16* {.importc: "__mmask16", imm.} = object
  mmask32* {.importc: "__mmask32", imm.} = object
  mmask64* {.importc: "__mmask64", imm.} = object

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
