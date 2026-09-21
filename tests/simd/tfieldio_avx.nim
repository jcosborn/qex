# tfieldio with the SSE,AVX composition: AVX512 removed from the build's setting.
# The pragmas must precede every import so the simd modules see the change.
when defined(AVX512):
  {.undef: AVX512.}
  include "../base/tfieldio"
else:
  echo "skipped: AVX512 is not in this build's simd setting"
