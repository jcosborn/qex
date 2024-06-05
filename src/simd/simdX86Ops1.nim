#import base
#getOptimPragmas()

proc baseImpl(b:NimNode; x:NimNode):NimNode =
  var n = x.len
  result = copyNimNode(x[0])
  template fold(b,x,y:untyped):untyped = b*x + y
  for i in 1..<n:
    result = getAst(fold(b, result, x[i]))
  #echo result.repr
#macro BASE(b:untyped; x:varargs[untyped]):auto = baseImpl(b, x)
macro BASE4(x:varargs[untyped]):auto = baseImpl(newLit(4), x)

# m128 operations

proc perm2*(x:m128):m128 {.alwaysInline.} =
  mm_shuffle_ps(x, x, BASE4(1,0,3,2).cuint)

proc perm1*(r:var m128; x:m128) {.alwaysInline.} =
  r = mm_shuffle_ps(x, x, BASE4(2,3,0,1).cuint)
proc perm2*(r:var m128; x:m128) {.alwaysInline.} =
  r = mm_shuffle_ps(x, x, BASE4(1,0,3,2).cuint)
proc perm4*(r:var m128; x:m128) {.alwaysInline.} =
  assert(false, "perm4 not valid for m128")
proc perm8*(r:var m128; x:m128) {.alwaysInline.} =
  assert(false, "perm8 not valid for m128")

var simdPermM128 = [
  mm_set_epi32(3,2,1,0),
  mm_set_epi32(2,3,0,1),
  mm_set_epi32(1,0,3,2),
  mm_set_epi32(0,1,2,3)
]
template perm*(x: m128, p: SomeNumber): untyped =
  mm_permutevar_ps(x, simdPermM128[p mod 4])

proc packp1*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  r[0] = t[1]
  l[1] = t[2]
  r[1] = t[3]
proc packm1*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  l[0] = t[1]
  r[1] = t[2]
  l[1] = t[3]
proc packp2*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  r[0] = t[2]
  r[1] = t[3]
proc packm2*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  l[0] = t[2]
  l[1] = t[3]
proc packp4*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packp4 not valid for m128")
proc packm4*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packm4 not valid for m128")
proc packp8*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packp8 not valid for m128")
proc packm8*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packm8 not valid for m128")

proc blendp1*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = r[0]
  t[2] = l[1]
  t[3] = r[1]
  assign(x, t)
proc blendm1*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = l[0]
  t[2] = r[1]
  t[3] = l[1]
  assign(x, t)
proc blendp2*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = l[1]
  t[2] = r[0]
  t[3] = r[1]
  assign(x, t)
proc blendm2*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = r[1]
  t[2] = l[0]
  t[3] = l[1]
  assign(x, t)
proc blendp4*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendp4 not valid for m128")
proc blendm4*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendm4 not valid for m128")
proc blendp8*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendp8 not valid for m128")
proc blendm8*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendm8 not valid for m128")


# m128d operations

proc perm1*(r: var m128d; x: m128d) {.alwaysInline.} =
  r = mm_shuffle_pd(x, x, cint(1))
proc perm2*(r: var m128d; x: m128d) {.alwaysInline.} =
  assert(false, "perm2 not valid for m128d")
proc perm4*(r: var m128d; x: m128d) {.alwaysInline.} =
  assert(false, "perm4 not valid for m128d")
proc perm8*(r: var m128d; x: m128d) {.alwaysInline.} =
  assert(false, "perm8 not valid for m128d")

proc packp1*(r: var openArray[SomeNumber]; x: m128d;
             l: var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  r[0] = t[1]
proc packm1*(r: var openArray[SomeNumber]; x: m128d;
             l: var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  l[0] = t[1]
proc packp2*(r: var openArray[SomeNumber]; x: m128d;
             l: var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packp2 not valid for m128d")
proc packm2*(r: var openArray[SomeNumber]; x: m128d;
             l: var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packp2 not valid for m128d")
proc packp4*(r: var openArray[SomeNumber]; x: m128d;
             l: var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packp4 not valid for m128d")
proc packm4*(r:var openArray[SomeNumber]; x:m128d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packm4 not valid for m128d")
proc packp8*(r:var openArray[SomeNumber]; x:m128d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packp8 not valid for m128d")
proc packm8*(r:var openArray[SomeNumber]; x:m128d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packm8 not valid for m128d")

proc blendp1*(x:var m128d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = r[0]
  assign(x, t)
proc blendm1*(x:var m128d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = l[0]
  assign(x, t)
proc blendp2*(x:var m128d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendp2 not valid for m128d")
proc blendm2*(x:var m128d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendm2 not valid for m128d")
proc blendp4*(x:var m128d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendp4 not valid for m128d")
proc blendm4*(x:var m128d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendm4 not valid for m128d")
proc blendp8*(x:var m128d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendp8 not valid for m128d")
proc blendm8*(x:var m128d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendm8 not valid for m128d")


# m256 operations

proc perm1*(r:var m256; x:m256) {.alwaysInline.} =
  r = mm256_permute_ps(x, BASE4(2,3,0,1))
proc perm2*(r:var m256; x:m256) {.alwaysInline.} =
  r = mm256_permute_ps(x, BASE4(1,0,3,2))
proc perm4*(r:var m256; x:m256) {.alwaysInline.} =
  r = mm256_permute2f128_ps(x, x, 1)
proc perm8*(r:var m256; x:m256) {.alwaysInline.} =
  assert(false, "perm8 not valid for m256")

var simdPermM256 = [
  mm256_set_epi32(7,6,5,4,3,2,1,0),
  mm256_set_epi32(6,7,4,5,2,3,0,1),
  mm256_set_epi32(5,4,7,6,1,0,3,2),
  mm256_set_epi32(4,5,6,7,0,1,2,3),
  mm256_set_epi32(3,2,1,0,7,6,5,4),
  mm256_set_epi32(2,3,0,1,6,7,4,5),
  mm256_set_epi32(1,0,3,2,5,4,7,6),
  mm256_set_epi32(0,1,2,3,4,5,6,7)
]
template perm*(x: m256, p: SomeNumber): untyped =
  mm256_permutevar8x32_ps(x, simdPermM256[p mod 8])

proc packp1*(r:var openArray[SomeNumber]; x:m256;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  r[0] = t[1]
  l[1] = t[2]
  r[1] = t[3]
  l[2] = t[4]
  r[2] = t[5]
  l[3] = t[6]
  r[3] = t[7]
proc packm1*(r:var openArray[SomeNumber]; x:m256;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  l[0] = t[1]
  r[1] = t[2]
  l[1] = t[3]
  r[2] = t[4]
  l[2] = t[5]
  r[3] = t[6]
  l[3] = t[7]
proc packp2*(r:var openArray[SomeNumber]; x:m256;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  r[0] = t[2]
  r[1] = t[3]
  l[2] = t[4]
  l[3] = t[5]
  r[2] = t[6]
  r[3] = t[7]
proc packm2*(r:var openArray[SomeNumber]; x:m256;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  l[0] = t[2]
  l[1] = t[3]
  r[2] = t[4]
  r[3] = t[5]
  l[2] = t[6]
  l[3] = t[7]
proc packp4*(r:var openArray[SomeNumber]; x:m256;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  l[2] = t[2]
  l[3] = t[3]
  r[0] = t[4]
  r[1] = t[5]
  r[2] = t[6]
  r[3] = t[7]
proc packm4*(r:var openArray[SomeNumber]; x:m256;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  r[2] = t[2]
  r[3] = t[3]
  l[0] = t[4]
  l[1] = t[5]
  l[2] = t[6]
  l[3] = t[7]
proc packp8*(r:var openArray[SomeNumber]; x:m256;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packp8 not valid for m256")
proc packm8*(r:var openArray[SomeNumber]; x:m256;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packm8 not valid for m256")

proc blendp1*(x:var m256; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = r[0]
  t[2] = l[1]
  t[3] = r[1]
  t[4] = l[2]
  t[5] = r[2]
  t[6] = l[3]
  t[7] = r[3]
  assign(x, t)
proc blendm1*(x:var m256; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = l[0]
  t[2] = r[1]
  t[3] = l[1]
  t[4] = r[2]
  t[5] = l[2]
  t[6] = r[3]
  t[7] = l[3]
  assign(x, t)
proc blendp2*(x:var m256; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = l[1]
  t[2] = r[0]
  t[3] = r[1]
  t[4] = l[2]
  t[5] = l[3]
  t[6] = r[2]
  t[7] = r[3]
  assign(x, t)
proc blendm2*(x:var m256; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = r[1]
  t[2] = l[0]
  t[3] = l[1]
  t[4] = r[2]
  t[5] = r[3]
  t[6] = l[2]
  t[7] = l[3]
  assign(x, t)
proc blendp4*(x:var m256; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = l[1]
  t[2] = l[2]
  t[3] = l[3]
  t[4] = r[0]
  t[5] = r[1]
  t[6] = r[2]
  t[7] = r[3]
  assign(x, t)
proc blendm4*(x:var m256; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = r[1]
  t[2] = r[2]
  t[3] = r[3]
  t[4] = l[0]
  t[5] = l[1]
  t[6] = l[2]
  t[7] = l[3]
  assign(x, t)
proc blendp8*(x:var m256; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendp8 not valid for m256")
proc blendm8*(x:var m256; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendm8 not valid for m256")


# m256d operations

var simdPermM256d = [
  mm256_set_epi32(7,6,5,4,3,2,1,0),
  mm256_set_epi32(5,4,7,6,1,0,3,2),
  mm256_set_epi32(3,2,1,0,7,6,5,4),
  mm256_set_epi32(1,0,3,2,5,4,7,6)
]
template perm*(x: m256d, p: SomeNumber): untyped =
  mm256_castps_pd(mm256_permutevar8x32_ps(mm256_castpd_ps(x),
                                          simdPermM256d[p mod 4]))

proc perm1*(r:var m256d; x:m256d) {.alwaysInline.} =
  r = mm256_permute_pd(x, 5)
proc perm2*(r:var m256d; x:m256d) {.alwaysInline.} =
  r = mm256_permute2f128_pd(x, x, 1)
proc perm4*(r:var m256d; x:m256d) {.alwaysInline.} =
  assert(false, "perm4 not valid for m256d")
proc perm8*(r:var m256d; x:m256d) {.alwaysInline.} =
  assert(false, "perm8 not valid for m256d")

proc packp1*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  r[0] = t[1]
  l[1] = t[2]
  r[1] = t[3]
proc packm1*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  l[0] = t[1]
  r[1] = t[2]
  l[1] = t[3]
proc packp2*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  r[0] = t[2]
  r[1] = t[3]
proc packm2*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  l[0] = t[2]
  l[1] = t[3]
proc packp4*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packp4 not valid for m256d")
proc packm4*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packm4 not valid for m256d")
proc packp8*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packp8 not valid for m256d")
proc packm8*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packm8 not valid for m256d")

proc blendp1*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = r[0]
  t[2] = l[1]
  t[3] = r[1]
  assign(x, t)
proc blendm1*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = l[0]
  t[2] = r[1]
  t[3] = l[1]
  assign(x, t)
proc blendp2*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = l[1]
  t[2] = r[0]
  t[3] = r[1]
  assign(x, t)
proc blendm2*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = r[1]
  t[2] = l[0]
  t[3] = l[1]
  assign(x, t)
proc blendp4*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendp4 not valid for m256d")
proc blendm4*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendm4 not valid for m256d")
proc blendp8*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendp8 not valid for m256d")
proc blendm8*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendm8 not valid for m256d")


# m512 operations

when defined(AVX512):
  var simdPermM512 = [
    mm512_set_epi32(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0)
    #mm512_set_epi32(6,7,4,5,2,3,0,1),
    #mm512_set_epi32(5,4,7,6,1,0,3,2),
    #mm512_set_epi32(4,5,6,7,0,1,2,3),
    #mm512_set_epi32(3,2,1,0,7,6,5,4),
    #mm512_set_epi32(2,3,0,1,6,7,4,5),
    #mm512_set_epi32(1,0,3,2,5,4,7,6),
    #mm512_set_epi32(0,1,2,3,4,5,6,7)
  ]
  template perm*(x: m512, p: SomeNumber): untyped =
    mm512_permutexvar_ps(simdPermM512[p mod 16], x)

#proc perm1*(r:var m512; x:m512) {.alwaysInline.} =
template perm1*(r:var m512; x:m512) =
  r = mm512_permute_ps(x, BASE4(2,3,0,1))
#proc perm2*(r:var m512; x:m512) {.alwaysInline.} =
template perm2*(r:var m512; x:m512) =
  r = mm512_permute_ps(x, BASE4(1,0,3,2))
#proc perm4*(r:var m512; x:m512) {.alwaysInline.} =
template perm4*(r:var m512; x:m512) =
  r = mm512_shuffle_f32x4(x, x, BASE4(2,3,0,1))
#proc perm8*(r:var m512; x:m512) {.alwaysInline.} =
template perm8*(r:var m512; x:m512) =
  r = mm512_shuffle_f32x4(x, x, BASE4(1,0,3,2))

proc packp1*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  r[0] = t[1]
  l[1] = t[2]
  r[1] = t[3]
  l[2] = t[4]
  r[2] = t[5]
  l[3] = t[6]
  r[3] = t[7]
  l[4] = t[8]
  r[4] = t[9]
  l[5] = t[10]
  r[5] = t[11]
  l[6] = t[12]
  r[6] = t[13]
  l[7] = t[14]
  r[7] = t[15]
proc packm1*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  l[0] = t[1]
  r[1] = t[2]
  l[1] = t[3]
  r[2] = t[4]
  l[2] = t[5]
  r[3] = t[6]
  l[3] = t[7]
  r[4] = t[8]
  l[4] = t[9]
  r[5] = t[10]
  l[5] = t[11]
  r[6] = t[12]
  l[6] = t[13]
  r[7] = t[14]
  l[7] = t[15]
proc packp2*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  r[0] = t[2]
  r[1] = t[3]
  l[2] = t[4]
  l[3] = t[5]
  r[2] = t[6]
  r[3] = t[7]
  l[4] = t[8]
  l[5] = t[9]
  r[4] = t[10]
  r[5] = t[11]
  l[6] = t[12]
  l[7] = t[13]
  r[6] = t[14]
  r[7] = t[15]
proc packm2*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  l[0] = t[2]
  l[1] = t[3]
  r[2] = t[4]
  r[3] = t[5]
  l[2] = t[6]
  l[3] = t[7]
  r[4] = t[8]
  r[5] = t[9]
  l[4] = t[10]
  l[5] = t[11]
  r[6] = t[12]
  r[7] = t[13]
  l[6] = t[14]
  l[7] = t[15]
proc packp4*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  l[2] = t[2]
  l[3] = t[3]
  r[0] = t[4]
  r[1] = t[5]
  r[2] = t[6]
  r[3] = t[7]
  l[4] = t[8]
  l[5] = t[9]
  l[6] = t[10]
  l[7] = t[11]
  r[4] = t[12]
  r[5] = t[13]
  r[6] = t[14]
  r[7] = t[15]
proc packm4*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  r[2] = t[2]
  r[3] = t[3]
  l[0] = t[4]
  l[1] = t[5]
  l[2] = t[6]
  l[3] = t[7]
  r[4] = t[8]
  r[5] = t[9]
  r[6] = t[10]
  r[7] = t[11]
  l[4] = t[12]
  l[5] = t[13]
  l[6] = t[14]
  l[7] = t[15]
proc packp8*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  l[2] = t[2]
  l[3] = t[3]
  l[4] = t[4]
  l[5] = t[5]
  l[6] = t[6]
  l[7] = t[7]
  r[0] = t[8]
  r[1] = t[9]
  r[2] = t[10]
  r[3] = t[11]
  r[4] = t[12]
  r[5] = t[13]
  r[6] = t[14]
  r[7] = t[15]
proc packm8*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  r[2] = t[2]
  r[3] = t[3]
  r[4] = t[4]
  r[5] = t[5]
  r[6] = t[6]
  r[7] = t[7]
  l[0] = t[8]
  l[1] = t[9]
  l[2] = t[10]
  l[3] = t[11]
  l[4] = t[12]
  l[5] = t[13]
  l[6] = t[14]
  l[7] = t[15]

proc blendp1*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0]  = l[0]
  t[1]  = r[0]
  t[2]  = l[1]
  t[3]  = r[1]
  t[4]  = l[2]
  t[5]  = r[2]
  t[6]  = l[3]
  t[7]  = r[3]
  t[8]  = l[4]
  t[9]  = r[4]
  t[10] = l[5]
  t[11] = r[5]
  t[12] = l[6]
  t[13] = r[6]
  t[14] = l[7]
  t[15] = r[7]
  assign(x, t)
proc blendm1*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0]  = r[0]
  t[1]  = l[0]
  t[2]  = r[1]
  t[3]  = l[1]
  t[4]  = r[2]
  t[5]  = l[2]
  t[6]  = r[3]
  t[7]  = l[3]
  t[8]  = r[4]
  t[9]  = l[4]
  t[10] = r[5]
  t[11] = l[5]
  t[12] = r[6]
  t[13] = l[6]
  t[14] = r[7]
  t[15] = l[7]
  assign(x, t)
proc blendp2*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0]  = l[0]
  t[1]  = l[1]
  t[2]  = r[0]
  t[3]  = r[1]
  t[4]  = l[2]
  t[5]  = l[3]
  t[6]  = r[2]
  t[7]  = r[3]
  t[8]  = l[4]
  t[9]  = l[5]
  t[10] = r[4]
  t[11] = r[5]
  t[12] = l[6]
  t[13] = l[7]
  t[14] = r[6]
  t[15] = r[7]
  assign(x, t)
proc blendm2*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0]  = r[0]
  t[1]  = r[1]
  t[2]  = l[0]
  t[3]  = l[1]
  t[4]  = r[2]
  t[5]  = r[3]
  t[6]  = l[2]
  t[7]  = l[3]
  t[8]  = r[4]
  t[9]  = r[5]
  t[10] = l[4]
  t[11] = l[5]
  t[12] = r[6]
  t[13] = r[7]
  t[14] = l[6]
  t[15] = l[7]
  assign(x, t)
proc blendp4*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0]  = l[0]
  t[1]  = l[1]
  t[2]  = l[2]
  t[3]  = l[3]
  t[4]  = r[0]
  t[5]  = r[1]
  t[6]  = r[2]
  t[7]  = r[3]
  t[8]  = l[4]
  t[9]  = l[5]
  t[10] = l[6]
  t[11] = l[7]
  t[12] = r[4]
  t[13] = r[5]
  t[14] = r[6]
  t[15] = r[7]
  assign(x, t)
proc blendm4*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0]  = r[0]
  t[1]  = r[1]
  t[2]  = r[2]
  t[3]  = r[3]
  t[4]  = l[0]
  t[5]  = l[1]
  t[6]  = l[2]
  t[7]  = l[3]
  t[8]  = r[4]
  t[9]  = r[5]
  t[10] = r[6]
  t[11] = r[7]
  t[12] = l[4]
  t[13] = l[5]
  t[14] = l[6]
  t[15] = l[7]
  assign(x, t)
proc blendp8*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0]  = l[0]
  t[1]  = l[1]
  t[2]  = l[2]
  t[3]  = l[3]
  t[4]  = l[4]
  t[5]  = l[5]
  t[6]  = l[6]
  t[7]  = l[7]
  t[8]  = r[0]
  t[9]  = r[1]
  t[10] = r[2]
  t[11] = r[3]
  t[12] = r[4]
  t[13] = r[5]
  t[14] = r[6]
  t[15] = r[7]
  assign(x, t)
proc blendm8*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0]  = r[0]
  t[1]  = r[1]
  t[2]  = r[2]
  t[3]  = r[3]
  t[4]  = r[4]
  t[5]  = r[5]
  t[6]  = r[6]
  t[7]  = r[7]
  t[8]  = l[0]
  t[9]  = l[1]
  t[10] = l[2]
  t[11] = l[3]
  t[12] = l[4]
  t[13] = l[5]
  t[14] = l[6]
  t[15] = l[7]
  assign(x, t)


# m512d operations

proc perm1*(r:var m512d; x:m512d) {.alwaysInline.} =
  r = mm512_permute_pd(x, BASE4(1,1,1,1))
proc perm2*(r:var m512d; x:m512d) {.alwaysInline.} =
  r = mm512_shuffle_f64x2(x, x, BASE4(2,3,0,1))
proc perm4*(r:var m512d; x:m512d) {.alwaysInline.} =
  r = mm512_shuffle_f64x2(x, x, BASE4(1,0,3,2))
proc perm8*(r:var m512d; x:m512d) {.alwaysInline.} =
  assert(false, "perm8 not valid for m512d")

proc packp1*(r:var openArray[SomeNumber]; x:m512d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  r[0] = t[1]
  l[1] = t[2]
  r[1] = t[3]
  l[2] = t[4]
  r[2] = t[5]
  l[3] = t[6]
  r[3] = t[7]
proc packm1*(r:var openArray[SomeNumber]; x:m512d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  l[0] = t[1]
  r[1] = t[2]
  l[1] = t[3]
  r[2] = t[4]
  l[2] = t[5]
  r[3] = t[6]
  l[3] = t[7]
proc packp2*(r:var openArray[SomeNumber]; x:m512d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  r[0] = t[2]
  r[1] = t[3]
  l[2] = t[4]
  l[3] = t[5]
  r[2] = t[6]
  r[3] = t[7]
proc packm2*(r:var openArray[SomeNumber]; x:m512d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  l[0] = t[2]
  l[1] = t[3]
  r[2] = t[4]
  r[3] = t[5]
  l[2] = t[6]
  l[3] = t[7]
proc packp4*(r:var openArray[SomeNumber]; x:m512d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  l[2] = t[2]
  l[3] = t[3]
  r[0] = t[4]
  r[1] = t[5]
  r[2] = t[6]
  r[3] = t[7]
proc packm4*(r:var openArray[SomeNumber]; x:m512d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  r[2] = t[2]
  r[3] = t[3]
  l[0] = t[4]
  l[1] = t[5]
  l[2] = t[6]
  l[3] = t[7]
proc packp8*(r:var openArray[SomeNumber]; x:m512d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packp8 not valid for m512d")
proc packm8*(r:var openArray[SomeNumber]; x:m512d;
             l:var openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "packm8 not valid for m512d")

proc blendp1*(x:var m512d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = r[0]
  t[2] = l[1]
  t[3] = r[1]
  t[4] = l[2]
  t[5] = r[2]
  t[6] = l[3]
  t[7] = r[3]
  assign(x, t)
proc blendm1*(x:var m512d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = l[0]
  t[2] = r[1]
  t[3] = l[1]
  t[4] = r[2]
  t[5] = l[2]
  t[6] = r[3]
  t[7] = l[3]
  assign(x, t)
proc blendp2*(x:var m512d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = l[1]
  t[2] = r[0]
  t[3] = r[1]
  t[4] = l[2]
  t[5] = l[3]
  t[6] = r[2]
  t[7] = r[3]
  assign(x, t)
proc blendm2*(x:var m512d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = r[1]
  t[2] = l[0]
  t[3] = l[1]
  t[4] = r[2]
  t[5] = r[3]
  t[6] = l[2]
  t[7] = l[3]
  assign(x, t)
proc blendp4*(x:var m512d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = l[1]
  t[2] = l[2]
  t[3] = l[3]
  t[4] = r[0]
  t[5] = r[1]
  t[6] = r[2]
  t[7] = r[3]
  assign(x, t)
proc blendm4*(x:var m512d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = r[1]
  t[2] = r[2]
  t[3] = r[3]
  t[4] = l[0]
  t[5] = l[1]
  t[6] = l[2]
  t[7] = l[3]
  assign(x, t)
proc blendp8*(x:var m512d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendp8 not valid for m512d")
proc blendm8*(x:var m512d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.alwaysInline.} =
  assert(false, "blendm8 not valid for m512d")
