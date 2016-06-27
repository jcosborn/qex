import macros
proc baseImpl(b:NimNode; x:NimNode):NimNode =
  var n = x.len
  result = copyNimNode(x[0])
  template fold(b,x,y:untyped):untyped = b*x + y
  for i in 1..<n:
    result = getAst(fold(b, result, x[i]))
  #echo result.repr
macro BASE(b:expr; x:varargs[expr]):auto = baseImpl(b, x)
macro BASE4(x:varargs[expr]):auto = baseImpl(newLit(4), x)

# m128 operations

proc perm1*(r:var m128; x:m128) {.inline.} =
  r = mm_shuffle_ps(x, x, BASE4(2,3,0,1).cuint)
proc perm2*(r:var m128; x:m128) {.inline.} =
  r = mm_shuffle_ps(x, x, BASE4(1,0,3,2).cuint)
proc perm4*(r:var m128; x:m128) {.inline.} =
  assert(false, "perm4 not valid for m128")

proc packp1*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  l[0] = t[0]
  r[0] = t[1]
  l[1] = t[2]
  r[1] = t[3]
proc packm1*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  l[0] = t[1]
  r[1] = t[2]
  l[1] = t[3]
proc packp2*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  r[0] = t[2]
  r[1] = t[3]
proc packm2*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  l[0] = t[2]
  l[1] = t[3]
proc packp4*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.inline.} =
  assert(false, "packp4 not valid for m128")
proc packm4*(r:var openArray[SomeNumber]; x:m128;
             l:var openArray[SomeNumber]) {.inline.} =
  assert(false, "packm4 not valid for m128")

proc blendp1*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = r[0]
  t[2] = l[1]
  t[3] = r[1]
  assign(x, t)
proc blendm1*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = l[0]
  t[2] = r[1]
  t[3] = l[1]
  assign(x, t)
proc blendp2*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = l[1]
  t[2] = r[0]
  t[3] = r[1]
  assign(x, t)
proc blendm2*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = r[1]
  t[2] = l[0]
  t[3] = l[1]
  assign(x, t)
proc blendp4*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  assert(false, "blendp4 not valid for m128")
proc blendm4*(x:var m128; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  assert(false, "blendm4 not valid for m128")

# m256 operations

proc perm1*(r:var m256; x:m256) {.inline.} =
  r = mm256_permute_ps(x, BASE4(2,3,0,1))
proc perm2*(r:var m256; x:m256) {.inline.} =
  r = mm256_permute_ps(x, BASE4(1,0,3,2))
proc perm4*(r:var m256; x:m256) {.inline.} =
  r = mm256_permute2f128_ps(x, x, 1)

proc packp1*(r:var openArray[SomeNumber]; x:m256;
             l:var openArray[SomeNumber]) {.inline.} =
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
             l:var openArray[SomeNumber]) {.inline.} =
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
             l:var openArray[SomeNumber]) {.inline.} =
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
             l:var openArray[SomeNumber]) {.inline.} =
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
             l:var openArray[SomeNumber]) {.inline.} =
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
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  r[2] = t[2]
  r[3] = t[3]
  l[0] = t[4]
  l[1] = t[5]
  l[2] = t[6]
  l[3] = t[7]

proc blendp1*(x:var m256; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
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
              l:openArray[SomeNumber]) {.inline.} =
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
              l:openArray[SomeNumber]) {.inline.} =
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
              l:openArray[SomeNumber]) {.inline.} =
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
              l:openArray[SomeNumber]) {.inline.} =
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
              l:openArray[SomeNumber]) {.inline.} =
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

# m256d operations

proc perm1*(r:var m256d; x:m256d) {.inline.} =
  r = mm256_permute_pd(x, 5)
proc perm2*(r:var m256d; x:m256d) {.inline.} =
  r = mm256_permute2f128_pd(x, x, 1)
proc perm4*(r:var m256d; x:m256d) {.inline.} =
  assert(false, "perm4 not valid for m256d")

proc packp1*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  l[0] = t[0]
  r[0] = t[1]
  l[1] = t[2]
  r[1] = t[3]
proc packm1*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  l[0] = t[1]
  r[1] = t[2]
  l[1] = t[3]
proc packp2*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  r[0] = t[2]
  r[1] = t[3]
proc packm2*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  l[0] = t[2]
  l[1] = t[3]
proc packp4*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.inline.} =
  assert(false, "packp4 not valid for m256d")
proc packm4*(r:var openArray[SomeNumber]; x:m256d;
             l:var openArray[SomeNumber]) {.inline.} =
  assert(false, "packm4 not valid for m256d")

proc blendp1*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = r[0]
  t[2] = l[1]
  t[3] = r[1]
  assign(x, t)
proc blendm1*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = l[0]
  t[2] = r[1]
  t[3] = l[1]
  assign(x, t)
proc blendp2*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = l[0]
  t[1] = l[1]
  t[2] = r[0]
  t[3] = r[1]
  assign(x, t)
proc blendm2*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  var t{.noInit.}:type(toArray(x))
  t[0] = r[0]
  t[1] = r[1]
  t[2] = l[0]
  t[3] = l[1]
  assign(x, t)
proc blendp4*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  assert(false, "blendp4 not valid for m256d")
proc blendm4*(x:var m256d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
  assert(false, "blendm4 not valid for m256d")


# m512 operations

proc perm1*(r:var m512; x:m512) {.inline.} =
  r = mm512_permute_ps(x, BASE4(2,3,0,1))
proc perm2*(r:var m512; x:m512) {.inline.} =
  r = mm512_permute_ps(x, BASE4(1,0,3,2))
proc perm4*(r:var m512; x:m512) {.inline.} =
  r = mm512_shuffle_f32x4(x, x, BASE4(2,3,0,1))
proc perm8*(r:var m512; x:m512) {.inline.} =
  r = mm512_shuffle_f32x4(x, x, BASE4(1,0,3,2))

proc packp1*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  l[0] = t[0]
  r[0] = t[1]
  l[1] = t[2]
  r[1] = t[3]
  l[2] = t[4]
  r[2] = t[5]
  l[3] = t[6]
  r[3] = t[7]
proc packm1*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  l[0] = t[1]
  r[1] = t[2]
  l[1] = t[3]
  r[2] = t[4]
  l[2] = t[5]
  r[3] = t[6]
  l[3] = t[7]
proc packp2*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  r[0] = t[2]
  r[1] = t[3]
  l[2] = t[4]
  l[3] = t[5]
  r[2] = t[6]
  r[3] = t[7]
proc packm2*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  l[0] = t[2]
  l[1] = t[3]
  r[2] = t[4]
  r[3] = t[5]
  l[2] = t[6]
  l[3] = t[7]
proc packp4*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  l[0] = t[0]
  l[1] = t[1]
  l[2] = t[2]
  l[3] = t[3]
  r[0] = t[4]
  r[1] = t[5]
  r[2] = t[6]
  r[3] = t[7]
proc packm4*(r:var openArray[SomeNumber]; x:m512;
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  r[2] = t[2]
  r[3] = t[3]
  l[0] = t[4]
  l[1] = t[5]
  l[2] = t[6]
  l[3] = t[7]

proc blendp1*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
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
proc blendm1*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
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
proc blendp2*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
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
proc blendm2*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
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
proc blendp4*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
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
proc blendm4*(x:var m512; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
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


# m512d operations

proc perm1*(r:var m512d; x:m512d) {.inline.} =
  r = mm512_permute_pd(x, BASE4(1,1,1,1))
proc perm2*(r:var m512d; x:m512d) {.inline.} =
  r = mm512_shuffle_f64x2(x, x, BASE4(2,3,0,1))
proc perm4*(r:var m512d; x:m512d) {.inline.} =
  r = mm512_shuffle_f64x2(x, x, BASE4(1,0,3,2))
proc perm8*(r:var m512d; x:m512d) {.inline.} =
  assert(false, "perm8 not valid for m512d")

proc packp1*(r:var openArray[SomeNumber]; x:m512d;
             l:var openArray[SomeNumber]) {.inline.} =
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
             l:var openArray[SomeNumber]) {.inline.} =
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
             l:var openArray[SomeNumber]) {.inline.} =
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
             l:var openArray[SomeNumber]) {.inline.} =
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
             l:var openArray[SomeNumber]) {.inline.} =
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
             l:var openArray[SomeNumber]) {.inline.} =
  let t = x.toArray
  r[0] = t[0]
  r[1] = t[1]
  r[2] = t[2]
  r[3] = t[3]
  l[0] = t[4]
  l[1] = t[5]
  l[2] = t[6]
  l[3] = t[7]

proc blendp1*(x:var m512d; r:openArray[SomeNumber];
              l:openArray[SomeNumber]) {.inline.} =
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
              l:openArray[SomeNumber]) {.inline.} =
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
              l:openArray[SomeNumber]) {.inline.} =
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
              l:openArray[SomeNumber]) {.inline.} =
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
              l:openArray[SomeNumber]) {.inline.} =
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
              l:openArray[SomeNumber]) {.inline.} =
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
