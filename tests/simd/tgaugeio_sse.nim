# tgaugeio with the SSE only composition: AVX and AVX512 removed from the build's setting.
# The pragmas must precede every import so the simd modules see the change.
when defined(AVX):
  {.undef: AVX512.}
  {.undef: AVX.}
  include "../base/tgaugeio"
else:
  echo "skipped: AVX is not in this build's simd setting"
