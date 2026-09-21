# tfieldio without intrinsics: every SIMD type becomes an array of scalars.
# The pragmas must precede every import so the simd modules see the change.
when defined(SSE):
  {.undef: AVX512.}
  {.undef: AVX.}
  {.undef: SSE.}
  include "../base/tfieldio"
else:
  echo "skipped: SSE is not in this build's simd setting"
