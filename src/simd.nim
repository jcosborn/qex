#import simdGcc
#export simdGcc

when defined(QPX):
  import simd/simdQpx
  export simdQpx
else:
  import simd/simdX86
  export simdX86
