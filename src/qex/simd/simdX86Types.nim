{.pragma: imm, header:"immintrin.h".}
type
  m64*   {.importc: "__m64"  , imm.} = object
  m128*  {.importc: "__m128" , imm.} = object
  m128d* {.importc: "__m128d", imm.} = object
  m128i* {.importc: "__m128i", imm.} = object
  m256*  {.importc: "__m256" , imm.} = object
  m256d* {.importc: "__m256d", imm.} = object
  m256i* {.importc: "__m256i", imm.} = object
  m512*  {.importc: "__m512" , imm.} = object
  m512d* {.importc: "__m512d", imm.} = object
  m512i* {.importc: "__m512i", imm.} = object
  mmask8*  {.importc: "__mmask8" , imm.} = object
  mmask16* {.importc: "__mmask16", imm.} = object
  mmask32* {.importc: "__mmask32", imm.} = object
  mmask64* {.importc: "__mmask64", imm.} = object

when defined(SSE):
  type
    SimdS4* = m128
    SimdD2* = m128d
    SimdI4* = m128i
    SimdH4* = distinct int64
when defined(AVX):
  type
    SimdS8* = m256
    SimdD4* = m256d
    SimdI8* = m256i
    SimdH8* = distinct SimdI4
when defined(AVX512):
  type
    SimdS16* = m512
    SimdD8*  = m512d
    SimdI16* = m512i
    SimdH16* = distinct SimdI8
