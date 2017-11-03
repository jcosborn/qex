##prefix _m_
##prefix _mm_
##prefix _MM_
##mangle _mm_load_ps load_m128
##mangle _mm_load_pd load_m128d
##mangle _mm_loadu_ps loadu_m128
##mangle _mm_loadu_pd loadu_m128d
##mangle _mm_load1_ps load1_m128
##mangle _mm_load1_pd load1_m128d
##mangle _mm_setzero_ps setzero_m128
##mangle _mm_setzero_pd setzero_m128d
##mangle _mm_set1_ps set1_m128
##mangle _mm_set1_pd set1_m128d

# this file was generated from simdSse.cnim with
# c2nim simdSse.cnim

import simdX86Types

proc mm_abs_epi16*(a: m128i): m128i {.importc: "_mm_abs_epi16", 
                                      header: "immintrin.h".}
proc mm_abs_epi32*(a: m128i): m128i {.importc: "_mm_abs_epi32", 
                                      header: "immintrin.h".}
proc mm_abs_epi8*(a: m128i): m128i {.importc: "_mm_abs_epi8", 
                                     header: "immintrin.h".}
proc mm_abs_pi16*(a: m64): m64 {.importc: "_mm_abs_pi16", header: "immintrin.h".}
proc mm_abs_pi32*(a: m64): m64 {.importc: "_mm_abs_pi32", header: "immintrin.h".}
proc mm_abs_pi8*(a: m64): m64 {.importc: "_mm_abs_pi8", header: "immintrin.h".}
proc mm_add_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_add_epi16", 
    header: "immintrin.h".}
proc mm_add_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_add_epi32", 
    header: "immintrin.h".}
proc mm_add_epi64*(a: m128i; b: m128i): m128i {.importc: "_mm_add_epi64", 
    header: "immintrin.h".}
proc mm_add_epi8*(a: m128i; b: m128i): m128i {.importc: "_mm_add_epi8", 
    header: "immintrin.h".}
proc mm_add_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_add_pd", 
    header: "immintrin.h".}
proc mm_add_pi16*(a: m64; b: m64): m64 {.importc: "_mm_add_pi16", 
    header: "immintrin.h".}
proc mm_add_pi32*(a: m64; b: m64): m64 {.importc: "_mm_add_pi32", 
    header: "immintrin.h".}
proc mm_add_pi8*(a: m64; b: m64): m64 {.importc: "_mm_add_pi8", 
                                        header: "immintrin.h".}
proc mm_add_ps*(a: m128; b: m128): m128 {.importc: "_mm_add_ps", 
    header: "immintrin.h".}
proc mm_add_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_add_sd", 
    header: "immintrin.h".}
proc mm_add_si64*(a: m64; b: m64): m64 {.importc: "_mm_add_si64", 
    header: "immintrin.h".}
proc mm_add_ss*(a: m128; b: m128): m128 {.importc: "_mm_add_ss", 
    header: "immintrin.h".}
proc mm_adds_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_adds_epi16", 
    header: "immintrin.h".}
proc mm_adds_epi8*(a: m128i; b: m128i): m128i {.importc: "_mm_adds_epi8", 
    header: "immintrin.h".}
proc mm_adds_epu16*(a: m128i; b: m128i): m128i {.importc: "_mm_adds_epu16", 
    header: "immintrin.h".}
proc mm_adds_epu8*(a: m128i; b: m128i): m128i {.importc: "_mm_adds_epu8", 
    header: "immintrin.h".}
proc mm_adds_pi16*(a: m64; b: m64): m64 {.importc: "_mm_adds_pi16", 
    header: "immintrin.h".}
proc mm_adds_pi8*(a: m64; b: m64): m64 {.importc: "_mm_adds_pi8", 
    header: "immintrin.h".}
proc mm_adds_pu16*(a: m64; b: m64): m64 {.importc: "_mm_adds_pu16", 
    header: "immintrin.h".}
proc mm_adds_pu8*(a: m64; b: m64): m64 {.importc: "_mm_adds_pu8", 
    header: "immintrin.h".}
proc mm_addsub_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_addsub_pd", 
    header: "immintrin.h".}
proc mm_addsub_ps*(a: m128; b: m128): m128 {.importc: "_mm_addsub_ps", 
    header: "immintrin.h".}
proc mm_alignr_epi8*(a: m128i; b: m128i; count: cint): m128i {.
    importc: "_mm_alignr_epi8", header: "immintrin.h".}
proc mm_alignr_pi8*(a: m64; b: m64; count: cint): m64 {.
    importc: "_mm_alignr_pi8", header: "immintrin.h".}
proc mm_and_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_and_pd", 
    header: "immintrin.h".}
proc mm_and_ps*(a: m128; b: m128): m128 {.importc: "_mm_and_ps", 
    header: "immintrin.h".}
proc mm_and_si128*(a: m128i; b: m128i): m128i {.importc: "_mm_and_si128", 
    header: "immintrin.h".}
proc mm_and_si64*(a: m64; b: m64): m64 {.importc: "_mm_and_si64", 
    header: "immintrin.h".}
proc mm_andnot_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_andnot_pd", 
    header: "immintrin.h".}
proc mm_andnot_ps*(a: m128; b: m128): m128 {.importc: "_mm_andnot_ps", 
    header: "immintrin.h".}
proc mm_andnot_si128*(a: m128i; b: m128i): m128i {.importc: "_mm_andnot_si128", 
    header: "immintrin.h".}
proc mm_andnot_si64*(a: m64; b: m64): m64 {.importc: "_mm_andnot_si64", 
    header: "immintrin.h".}
proc mm_avg_epu16*(a: m128i; b: m128i): m128i {.importc: "_mm_avg_epu16", 
    header: "immintrin.h".}
proc mm_avg_epu8*(a: m128i; b: m128i): m128i {.importc: "_mm_avg_epu8", 
    header: "immintrin.h".}
proc mm_avg_pu16*(a: m64; b: m64): m64 {.importc: "_mm_avg_pu16", 
    header: "immintrin.h".}
proc mm_avg_pu8*(a: m64; b: m64): m64 {.importc: "_mm_avg_pu8", 
                                        header: "immintrin.h".}
proc mm_blend_epi16*(a: m128i; b: m128i; imm8: cint): m128i {.
    importc: "_mm_blend_epi16", header: "immintrin.h".}
proc mm_blend_pd*(a: m128d; b: m128d; imm8: cint): m128d {.
    importc: "_mm_blend_pd", header: "immintrin.h".}
proc mm_blend_ps*(a: m128; b: m128; imm8: cint): m128 {.importc: "_mm_blend_ps", 
    header: "immintrin.h".}
proc mm_blendv_epi8*(a: m128i; b: m128i; mask: m128i): m128i {.
    importc: "_mm_blendv_epi8", header: "immintrin.h".}
proc mm_blendv_pd*(a: m128d; b: m128d; mask: m128d): m128d {.
    importc: "_mm_blendv_pd", header: "immintrin.h".}
proc mm_blendv_ps*(a: m128; b: m128; mask: m128): m128 {.
    importc: "_mm_blendv_ps", header: "immintrin.h".}
proc mm_bslli_si128*(a: m128i; imm8: cint): m128i {.importc: "_mm_bslli_si128", 
    header: "immintrin.h".}
proc mm_bsrli_si128*(a: m128i; imm8: cint): m128i {.importc: "_mm_bsrli_si128", 
    header: "immintrin.h".}
proc mm_castpd_ps*(a: m128d): m128 {.importc: "_mm_castpd_ps", 
                                     header: "immintrin.h".}
proc mm_castpd_si128*(a: m128d): m128i {.importc: "_mm_castpd_si128", 
    header: "immintrin.h".}
proc mm_castps_pd*(a: m128): m128d {.importc: "_mm_castps_pd", 
                                     header: "immintrin.h".}
proc mm_castps_si128*(a: m128): m128i {.importc: "_mm_castps_si128", 
                                        header: "immintrin.h".}
proc mm_castsi128_pd*(a: m128i): m128d {.importc: "_mm_castsi128_pd", 
    header: "immintrin.h".}
proc mm_castsi128_ps*(a: m128i): m128 {.importc: "_mm_castsi128_ps", 
                                        header: "immintrin.h".}
proc mm_ceil_pd*(a: m128d): m128d {.importc: "_mm_ceil_pd", 
                                    header: "immintrin.h".}
proc mm_ceil_ps*(a: m128): m128 {.importc: "_mm_ceil_ps", header: "immintrin.h".}
proc mm_ceil_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_ceil_sd", 
    header: "immintrin.h".}
proc mm_ceil_ss*(a: m128; b: m128): m128 {.importc: "_mm_ceil_ss", 
    header: "immintrin.h".}
proc mm_clflush*(p: pointer) {.importc: "_mm_clflush", header: "immintrin.h".}
proc mm_cmpeq_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_cmpeq_epi16", 
    header: "immintrin.h".}
proc mm_cmpeq_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_cmpeq_epi32", 
    header: "immintrin.h".}
proc mm_cmpeq_epi64*(a: m128i; b: m128i): m128i {.importc: "_mm_cmpeq_epi64", 
    header: "immintrin.h".}
proc mm_cmpeq_epi8*(a: m128i; b: m128i): m128i {.importc: "_mm_cmpeq_epi8", 
    header: "immintrin.h".}
proc mm_cmpeq_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpeq_pd", 
    header: "immintrin.h".}
proc mm_cmpeq_pi16*(a: m64; b: m64): m64 {.importc: "_mm_cmpeq_pi16", 
    header: "immintrin.h".}
proc mm_cmpeq_pi32*(a: m64; b: m64): m64 {.importc: "_mm_cmpeq_pi32", 
    header: "immintrin.h".}
proc mm_cmpeq_pi8*(a: m64; b: m64): m64 {.importc: "_mm_cmpeq_pi8", 
    header: "immintrin.h".}
proc mm_cmpeq_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmpeq_ps", 
    header: "immintrin.h".}
proc mm_cmpeq_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpeq_sd", 
    header: "immintrin.h".}
proc mm_cmpeq_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmpeq_ss", 
    header: "immintrin.h".}
proc mm_cmpestra*(a: m128i; la: cint; b: m128i; lb: cint; imm8: cint): cint {.
    importc: "_mm_cmpestra", header: "immintrin.h".}
proc mm_cmpestrc*(a: m128i; la: cint; b: m128i; lb: cint; imm8: cint): cint {.
    importc: "_mm_cmpestrc", header: "immintrin.h".}
proc mm_cmpestri*(a: m128i; la: cint; b: m128i; lb: cint; imm8: cint): cint {.
    importc: "_mm_cmpestri", header: "immintrin.h".}
proc mm_cmpestrm*(a: m128i; la: cint; b: m128i; lb: cint; imm8: cint): m128i {.
    importc: "_mm_cmpestrm", header: "immintrin.h".}
proc mm_cmpestro*(a: m128i; la: cint; b: m128i; lb: cint; imm8: cint): cint {.
    importc: "_mm_cmpestro", header: "immintrin.h".}
proc mm_cmpestrs*(a: m128i; la: cint; b: m128i; lb: cint; imm8: cint): cint {.
    importc: "_mm_cmpestrs", header: "immintrin.h".}
proc mm_cmpestrz*(a: m128i; la: cint; b: m128i; lb: cint; imm8: cint): cint {.
    importc: "_mm_cmpestrz", header: "immintrin.h".}
proc mm_cmpge_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpge_pd", 
    header: "immintrin.h".}
proc mm_cmpge_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmpge_ps", 
    header: "immintrin.h".}
proc mm_cmpge_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpge_sd", 
    header: "immintrin.h".}
proc mm_cmpge_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmpge_ss", 
    header: "immintrin.h".}
proc mm_cmpgt_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_cmpgt_epi16", 
    header: "immintrin.h".}
proc mm_cmpgt_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_cmpgt_epi32", 
    header: "immintrin.h".}
proc mm_cmpgt_epi64*(a: m128i; b: m128i): m128i {.importc: "_mm_cmpgt_epi64", 
    header: "immintrin.h".}
proc mm_cmpgt_epi8*(a: m128i; b: m128i): m128i {.importc: "_mm_cmpgt_epi8", 
    header: "immintrin.h".}
proc mm_cmpgt_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpgt_pd", 
    header: "immintrin.h".}
proc mm_cmpgt_pi16*(a: m64; b: m64): m64 {.importc: "_mm_cmpgt_pi16", 
    header: "immintrin.h".}
proc mm_cmpgt_pi32*(a: m64; b: m64): m64 {.importc: "_mm_cmpgt_pi32", 
    header: "immintrin.h".}
proc mm_cmpgt_pi8*(a: m64; b: m64): m64 {.importc: "_mm_cmpgt_pi8", 
    header: "immintrin.h".}
proc mm_cmpgt_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmpgt_ps", 
    header: "immintrin.h".}
proc mm_cmpgt_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpgt_sd", 
    header: "immintrin.h".}
proc mm_cmpgt_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmpgt_ss", 
    header: "immintrin.h".}
proc mm_cmpistra*(a: m128i; b: m128i; imm8: cint): cint {.
    importc: "_mm_cmpistra", header: "immintrin.h".}
proc mm_cmpistrc*(a: m128i; b: m128i; imm8: cint): cint {.
    importc: "_mm_cmpistrc", header: "immintrin.h".}
proc mm_cmpistri*(a: m128i; b: m128i; imm8: cint): cint {.
    importc: "_mm_cmpistri", header: "immintrin.h".}
proc mm_cmpistrm*(a: m128i; b: m128i; imm8: cint): m128i {.
    importc: "_mm_cmpistrm", header: "immintrin.h".}
proc mm_cmpistro*(a: m128i; b: m128i; imm8: cint): cint {.
    importc: "_mm_cmpistro", header: "immintrin.h".}
proc mm_cmpistrs*(a: m128i; b: m128i; imm8: cint): cint {.
    importc: "_mm_cmpistrs", header: "immintrin.h".}
proc mm_cmpistrz*(a: m128i; b: m128i; imm8: cint): cint {.
    importc: "_mm_cmpistrz", header: "immintrin.h".}
proc mm_cmple_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmple_pd", 
    header: "immintrin.h".}
proc mm_cmple_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmple_ps", 
    header: "immintrin.h".}
proc mm_cmple_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmple_sd", 
    header: "immintrin.h".}
proc mm_cmple_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmple_ss", 
    header: "immintrin.h".}
proc mm_cmplt_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_cmplt_epi16", 
    header: "immintrin.h".}
proc mm_cmplt_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_cmplt_epi32", 
    header: "immintrin.h".}
proc mm_cmplt_epi8*(a: m128i; b: m128i): m128i {.importc: "_mm_cmplt_epi8", 
    header: "immintrin.h".}
proc mm_cmplt_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmplt_pd", 
    header: "immintrin.h".}
proc mm_cmplt_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmplt_ps", 
    header: "immintrin.h".}
proc mm_cmplt_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmplt_sd", 
    header: "immintrin.h".}
proc mm_cmplt_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmplt_ss", 
    header: "immintrin.h".}
proc mm_cmpneq_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpneq_pd", 
    header: "immintrin.h".}
proc mm_cmpneq_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmpneq_ps", 
    header: "immintrin.h".}
proc mm_cmpneq_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpneq_sd", 
    header: "immintrin.h".}
proc mm_cmpneq_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmpneq_ss", 
    header: "immintrin.h".}
proc mm_cmpnge_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpnge_pd", 
    header: "immintrin.h".}
proc mm_cmpnge_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmpnge_ps", 
    header: "immintrin.h".}
proc mm_cmpnge_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpnge_sd", 
    header: "immintrin.h".}
proc mm_cmpnge_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmpnge_ss", 
    header: "immintrin.h".}
proc mm_cmpngt_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpngt_pd", 
    header: "immintrin.h".}
proc mm_cmpngt_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmpngt_ps", 
    header: "immintrin.h".}
proc mm_cmpngt_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpngt_sd", 
    header: "immintrin.h".}
proc mm_cmpngt_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmpngt_ss", 
    header: "immintrin.h".}
proc mm_cmpnle_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpnle_pd", 
    header: "immintrin.h".}
proc mm_cmpnle_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmpnle_ps", 
    header: "immintrin.h".}
proc mm_cmpnle_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpnle_sd", 
    header: "immintrin.h".}
proc mm_cmpnle_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmpnle_ss", 
    header: "immintrin.h".}
proc mm_cmpnlt_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpnlt_pd", 
    header: "immintrin.h".}
proc mm_cmpnlt_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmpnlt_ps", 
    header: "immintrin.h".}
proc mm_cmpnlt_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpnlt_sd", 
    header: "immintrin.h".}
proc mm_cmpnlt_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmpnlt_ss", 
    header: "immintrin.h".}
proc mm_cmpord_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpord_pd", 
    header: "immintrin.h".}
proc mm_cmpord_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmpord_ps", 
    header: "immintrin.h".}
proc mm_cmpord_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpord_sd", 
    header: "immintrin.h".}
proc mm_cmpord_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmpord_ss", 
    header: "immintrin.h".}
proc mm_cmpunord_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpunord_pd", 
    header: "immintrin.h".}
proc mm_cmpunord_ps*(a: m128; b: m128): m128 {.importc: "_mm_cmpunord_ps", 
    header: "immintrin.h".}
proc mm_cmpunord_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_cmpunord_sd", 
    header: "immintrin.h".}
proc mm_cmpunord_ss*(a: m128; b: m128): m128 {.importc: "_mm_cmpunord_ss", 
    header: "immintrin.h".}
proc mm_comieq_sd*(a: m128d; b: m128d): cint {.importc: "_mm_comieq_sd", 
    header: "immintrin.h".}
proc mm_comieq_ss*(a: m128; b: m128): cint {.importc: "_mm_comieq_ss", 
    header: "immintrin.h".}
proc mm_comige_sd*(a: m128d; b: m128d): cint {.importc: "_mm_comige_sd", 
    header: "immintrin.h".}
proc mm_comige_ss*(a: m128; b: m128): cint {.importc: "_mm_comige_ss", 
    header: "immintrin.h".}
proc mm_comigt_sd*(a: m128d; b: m128d): cint {.importc: "_mm_comigt_sd", 
    header: "immintrin.h".}
proc mm_comigt_ss*(a: m128; b: m128): cint {.importc: "_mm_comigt_ss", 
    header: "immintrin.h".}
proc mm_comile_sd*(a: m128d; b: m128d): cint {.importc: "_mm_comile_sd", 
    header: "immintrin.h".}
proc mm_comile_ss*(a: m128; b: m128): cint {.importc: "_mm_comile_ss", 
    header: "immintrin.h".}
proc mm_comilt_sd*(a: m128d; b: m128d): cint {.importc: "_mm_comilt_sd", 
    header: "immintrin.h".}
proc mm_comilt_ss*(a: m128; b: m128): cint {.importc: "_mm_comilt_ss", 
    header: "immintrin.h".}
proc mm_comineq_sd*(a: m128d; b: m128d): cint {.importc: "_mm_comineq_sd", 
    header: "immintrin.h".}
proc mm_comineq_ss*(a: m128; b: m128): cint {.importc: "_mm_comineq_ss", 
    header: "immintrin.h".}
proc mm_crc32_u16*(crc: cuint; v: cushort): cuint {.importc: "_mm_crc32_u16", 
    header: "immintrin.h".}
proc mm_crc32_u32*(crc: cuint; v: cuint): cuint {.importc: "_mm_crc32_u32", 
    header: "immintrin.h".}
proc mm_crc32_u64*(crc: uint64; v: uint64): uint64 {.importc: "_mm_crc32_u64", 
    header: "immintrin.h".}
proc mm_crc32_u8*(crc: cuint; v: cuchar): cuint {.importc: "_mm_crc32_u8", 
    header: "immintrin.h".}
proc mm_cvt_pi2ps*(a: m128; b: m64): m128 {.importc: "_mm_cvt_pi2ps", 
    header: "immintrin.h".}
proc mm_cvt_ps2pi*(a: m128): m64 {.importc: "_mm_cvt_ps2pi", 
                                   header: "immintrin.h".}
proc mm_cvt_si2ss*(a: m128; b: cint): m128 {.importc: "_mm_cvt_si2ss", 
    header: "immintrin.h".}
proc mm_cvt_ss2si*(a: m128): cint {.importc: "_mm_cvt_ss2si", 
                                    header: "immintrin.h".}
proc mm_cvtepi16_epi32*(a: m128i): m128i {.importc: "_mm_cvtepi16_epi32", 
    header: "immintrin.h".}
proc mm_cvtepi16_epi64*(a: m128i): m128i {.importc: "_mm_cvtepi16_epi64", 
    header: "immintrin.h".}
proc mm_cvtepi32_epi64*(a: m128i): m128i {.importc: "_mm_cvtepi32_epi64", 
    header: "immintrin.h".}
proc mm_cvtepi32_pd*(a: m128i): m128d {.importc: "_mm_cvtepi32_pd", 
                                        header: "immintrin.h".}
proc mm_cvtepi32_ps*(a: m128i): m128 {.importc: "_mm_cvtepi32_ps", 
                                       header: "immintrin.h".}
proc mm_cvtepi8_epi16*(a: m128i): m128i {.importc: "_mm_cvtepi8_epi16", 
    header: "immintrin.h".}
proc mm_cvtepi8_epi32*(a: m128i): m128i {.importc: "_mm_cvtepi8_epi32", 
    header: "immintrin.h".}
proc mm_cvtepi8_epi64*(a: m128i): m128i {.importc: "_mm_cvtepi8_epi64", 
    header: "immintrin.h".}
proc mm_cvtepu16_epi32*(a: m128i): m128i {.importc: "_mm_cvtepu16_epi32", 
    header: "immintrin.h".}
proc mm_cvtepu16_epi64*(a: m128i): m128i {.importc: "_mm_cvtepu16_epi64", 
    header: "immintrin.h".}
proc mm_cvtepu32_epi64*(a: m128i): m128i {.importc: "_mm_cvtepu32_epi64", 
    header: "immintrin.h".}
proc mm_cvtepu8_epi16*(a: m128i): m128i {.importc: "_mm_cvtepu8_epi16", 
    header: "immintrin.h".}
proc mm_cvtepu8_epi32*(a: m128i): m128i {.importc: "_mm_cvtepu8_epi32", 
    header: "immintrin.h".}
proc mm_cvtepu8_epi64*(a: m128i): m128i {.importc: "_mm_cvtepu8_epi64", 
    header: "immintrin.h".}
proc mm_cvtm64_si64*(a: m64): int64 {.importc: "_mm_cvtm64_si64", 
                                      header: "immintrin.h".}
proc mm_cvtpd_epi32*(a: m128d): m128i {.importc: "_mm_cvtpd_epi32", 
                                        header: "immintrin.h".}
proc mm_cvtpd_pi32*(a: m128d): m64 {.importc: "_mm_cvtpd_pi32", 
                                     header: "immintrin.h".}
proc mm_cvtpd_ps*(a: m128d): m128 {.importc: "_mm_cvtpd_ps", 
                                    header: "immintrin.h".}
proc mm_cvtpi16_ps*(a: m64): m128 {.importc: "_mm_cvtpi16_ps", 
                                    header: "immintrin.h".}
proc mm_cvtpi32_pd*(a: m64): m128d {.importc: "_mm_cvtpi32_pd", 
                                     header: "immintrin.h".}
proc mm_cvtpi32_ps*(a: m128; b: m64): m128 {.importc: "_mm_cvtpi32_ps", 
    header: "immintrin.h".}
proc mm_cvtpi32x2_ps*(a: m64; b: m64): m128 {.importc: "_mm_cvtpi32x2_ps", 
    header: "immintrin.h".}
proc mm_cvtpi8_ps*(a: m64): m128 {.importc: "_mm_cvtpi8_ps", 
                                   header: "immintrin.h".}
proc mm_cvtps_epi32*(a: m128): m128i {.importc: "_mm_cvtps_epi32", 
                                       header: "immintrin.h".}
proc mm_cvtps_pd*(a: m128): m128d {.importc: "_mm_cvtps_pd", 
                                    header: "immintrin.h".}
proc mm_cvtps_pi16*(a: m128): m64 {.importc: "_mm_cvtps_pi16", 
                                    header: "immintrin.h".}
proc mm_cvtps_pi32*(a: m128): m64 {.importc: "_mm_cvtps_pi32", 
                                    header: "immintrin.h".}
proc mm_cvtps_pi8*(a: m128): m64 {.importc: "_mm_cvtps_pi8", 
                                   header: "immintrin.h".}
proc mm_cvtpu16_ps*(a: m64): m128 {.importc: "_mm_cvtpu16_ps", 
                                    header: "immintrin.h".}
proc mm_cvtpu8_ps*(a: m64): m128 {.importc: "_mm_cvtpu8_ps", 
                                   header: "immintrin.h".}
proc mm_cvtsd_f64*(a: m128d): cdouble {.importc: "_mm_cvtsd_f64", 
                                        header: "immintrin.h".}
proc mm_cvtsd_si32*(a: m128d): cint {.importc: "_mm_cvtsd_si32", 
                                      header: "immintrin.h".}
proc mm_cvtsd_si64*(a: m128d): int64 {.importc: "_mm_cvtsd_si64", 
                                       header: "immintrin.h".}
proc mm_cvtsd_si64x*(a: m128d): int64 {.importc: "_mm_cvtsd_si64x", 
                                        header: "immintrin.h".}
proc mm_cvtsd_ss*(a: m128; b: m128d): m128 {.importc: "_mm_cvtsd_ss", 
    header: "immintrin.h".}
proc mm_cvtsi128_si32*(a: m128i): cint {.importc: "_mm_cvtsi128_si32", 
    header: "immintrin.h".}
proc mm_cvtsi128_si64*(a: m128i): int64 {.importc: "_mm_cvtsi128_si64", 
    header: "immintrin.h".}
proc mm_cvtsi128_si64x*(a: m128i): int64 {.importc: "_mm_cvtsi128_si64x", 
    header: "immintrin.h".}
proc mm_cvtsi32_sd*(a: m128d; b: cint): m128d {.importc: "_mm_cvtsi32_sd", 
    header: "immintrin.h".}
proc mm_cvtsi32_si128*(a: cint): m128i {.importc: "_mm_cvtsi32_si128", 
    header: "immintrin.h".}
proc mm_cvtsi32_si64*(a: cint): m64 {.importc: "_mm_cvtsi32_si64", 
                                      header: "immintrin.h".}
proc mm_cvtsi32_ss*(a: m128; b: cint): m128 {.importc: "_mm_cvtsi32_ss", 
    header: "immintrin.h".}
proc mm_cvtsi64_m64*(a: int64): m64 {.importc: "_mm_cvtsi64_m64", 
                                      header: "immintrin.h".}
proc mm_cvtsi64_sd*(a: m128d; b: int64): m128d {.importc: "_mm_cvtsi64_sd", 
    header: "immintrin.h".}
proc mm_cvtsi64_si128*(a: int64): m128i {.importc: "_mm_cvtsi64_si128", 
    header: "immintrin.h".}
proc mm_cvtsi64_si32*(a: m64): cint {.importc: "_mm_cvtsi64_si32", 
                                      header: "immintrin.h".}
proc mm_cvtsi64_ss*(a: m128; b: int64): m128 {.importc: "_mm_cvtsi64_ss", 
    header: "immintrin.h".}
proc mm_cvtsi64x_sd*(a: m128d; b: int64): m128d {.importc: "_mm_cvtsi64x_sd", 
    header: "immintrin.h".}
proc mm_cvtsi64x_si128*(a: int64): m128i {.importc: "_mm_cvtsi64x_si128", 
    header: "immintrin.h".}
proc mm_cvtss_f32*(a: m128): cfloat {.importc: "_mm_cvtss_f32", 
                                      header: "immintrin.h".}
proc mm_cvtss_sd*(a: m128d; b: m128): m128d {.importc: "_mm_cvtss_sd", 
    header: "immintrin.h".}
proc mm_cvtss_si32*(a: m128): cint {.importc: "_mm_cvtss_si32", 
                                     header: "immintrin.h".}
proc mm_cvtss_si64*(a: m128): int64 {.importc: "_mm_cvtss_si64", 
                                      header: "immintrin.h".}
proc mm_cvtt_ps2pi*(a: m128): m64 {.importc: "_mm_cvtt_ps2pi", 
                                    header: "immintrin.h".}
proc mm_cvtt_ss2si*(a: m128): cint {.importc: "_mm_cvtt_ss2si", 
                                     header: "immintrin.h".}
proc mm_cvttpd_epi32*(a: m128d): m128i {.importc: "_mm_cvttpd_epi32", 
    header: "immintrin.h".}
proc mm_cvttpd_pi32*(a: m128d): m64 {.importc: "_mm_cvttpd_pi32", 
                                      header: "immintrin.h".}
proc mm_cvttps_epi32*(a: m128): m128i {.importc: "_mm_cvttps_epi32", 
                                        header: "immintrin.h".}
proc mm_cvttps_pi32*(a: m128): m64 {.importc: "_mm_cvttps_pi32", 
                                     header: "immintrin.h".}
proc mm_cvttsd_si32*(a: m128d): cint {.importc: "_mm_cvttsd_si32", 
                                       header: "immintrin.h".}
proc mm_cvttsd_si64*(a: m128d): int64 {.importc: "_mm_cvttsd_si64", 
                                        header: "immintrin.h".}
proc mm_cvttsd_si64x*(a: m128d): int64 {.importc: "_mm_cvttsd_si64x", 
    header: "immintrin.h".}
proc mm_cvttss_si32*(a: m128): cint {.importc: "_mm_cvttss_si32", 
                                      header: "immintrin.h".}
proc mm_cvttss_si64*(a: m128): int64 {.importc: "_mm_cvttss_si64", 
                                       header: "immintrin.h".}
proc mm_div_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_div_pd", 
    header: "immintrin.h".}
proc mm_div_ps*(a: m128; b: m128): m128 {.importc: "_mm_div_ps", 
    header: "immintrin.h".}
proc mm_div_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_div_sd", 
    header: "immintrin.h".}
proc mm_div_ss*(a: m128; b: m128): m128 {.importc: "_mm_div_ss", 
    header: "immintrin.h".}
proc mm_dp_pd*(a: m128d; b: m128d; imm8: cint): m128d {.importc: "_mm_dp_pd", 
    header: "immintrin.h".}
proc mm_dp_ps*(a: m128; b: m128; imm8: cint): m128 {.importc: "_mm_dp_ps", 
    header: "immintrin.h".}
#redundant: void _m_empty (void);

proc mm_empty*() {.importc: "_mm_empty", header: "immintrin.h".}
proc mm_extract_epi16*(a: m128i; imm8: cint): cint {.
    importc: "_mm_extract_epi16", header: "immintrin.h".}
proc mm_extract_epi32*(a: m128i; imm8: cint): cint {.
    importc: "_mm_extract_epi32", header: "immintrin.h".}
proc mm_extract_epi64*(a: m128i; imm8: cint): int64 {.
    importc: "_mm_extract_epi64", header: "immintrin.h".}
proc mm_extract_epi8*(a: m128i; imm8: cint): cint {.importc: "_mm_extract_epi8", 
    header: "immintrin.h".}
proc mm_extract_pi16*(a: m64; imm8: cint): cint {.importc: "_mm_extract_pi16", 
    header: "immintrin.h".}
proc mm_extract_ps*(a: m128; imm8: cint): cint {.importc: "_mm_extract_ps", 
    header: "immintrin.h".}
proc mm_floor_pd*(a: m128d): m128d {.importc: "_mm_floor_pd", 
                                     header: "immintrin.h".}
proc mm_floor_ps*(a: m128): m128 {.importc: "_mm_floor_ps", 
                                   header: "immintrin.h".}
proc mm_floor_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_floor_sd", 
    header: "immintrin.h".}
proc mm_floor_ss*(a: m128; b: m128): m128 {.importc: "_mm_floor_ss", 
    header: "immintrin.h".}
proc m_from_int*(a: cint): m64 {.importc: "_m_from_int", header: "immintrin.h".}
proc m_from_int64*(a: int64): m64 {.importc: "_m_from_int64", 
                                    header: "immintrin.h".}
proc MM_GET_EXCEPTION_MASK*(): cuint {.importc: "_MM_GET_EXCEPTION_MASK", 
                                       header: "immintrin.h".}
proc MM_GET_EXCEPTION_STATE*(): cuint {.importc: "_MM_GET_EXCEPTION_STATE", 
                                        header: "immintrin.h".}
proc MM_GET_FLUSH_ZERO_MODE*(): cuint {.importc: "_MM_GET_FLUSH_ZERO_MODE", 
                                        header: "immintrin.h".}
proc MM_GET_ROUNDING_MODE*(): cuint {.importc: "_MM_GET_ROUNDING_MODE", 
                                      header: "immintrin.h".}
proc mm_getcsr*(): cuint {.importc: "_mm_getcsr", header: "immintrin.h".}
proc mm_hadd_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_hadd_epi16", 
    header: "immintrin.h".}
proc mm_hadd_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_hadd_epi32", 
    header: "immintrin.h".}
proc mm_hadd_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_hadd_pd", 
    header: "immintrin.h".}
proc mm_hadd_pi16*(a: m64; b: m64): m64 {.importc: "_mm_hadd_pi16", 
    header: "immintrin.h".}
proc mm_hadd_pi32*(a: m64; b: m64): m64 {.importc: "_mm_hadd_pi32", 
    header: "immintrin.h".}
proc mm_hadd_ps*(a: m128; b: m128): m128 {.importc: "_mm_hadd_ps", 
    header: "immintrin.h".}
proc mm_hadds_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_hadds_epi16", 
    header: "immintrin.h".}
proc mm_hadds_pi16*(a: m64; b: m64): m64 {.importc: "_mm_hadds_pi16", 
    header: "immintrin.h".}
proc mm_hsub_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_hsub_epi16", 
    header: "immintrin.h".}
proc mm_hsub_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_hsub_epi32", 
    header: "immintrin.h".}
proc mm_hsub_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_hsub_pd", 
    header: "immintrin.h".}
proc mm_hsub_pi16*(a: m64; b: m64): m64 {.importc: "_mm_hsub_pi16", 
    header: "immintrin.h".}
proc mm_hsub_pi32*(a: m64; b: m64): m64 {.importc: "_mm_hsub_pi32", 
    header: "immintrin.h".}
proc mm_hsub_ps*(a: m128; b: m128): m128 {.importc: "_mm_hsub_ps", 
    header: "immintrin.h".}
proc mm_hsubs_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_hsubs_epi16", 
    header: "immintrin.h".}
proc mm_hsubs_pi16*(a: m64; b: m64): m64 {.importc: "_mm_hsubs_pi16", 
    header: "immintrin.h".}
proc mm_insert_epi16*(a: m128i; i: cint; imm8: cint): m128i {.
    importc: "_mm_insert_epi16", header: "immintrin.h".}
proc mm_insert_epi32*(a: m128i; i: cint; imm8: cint): m128i {.
    importc: "_mm_insert_epi32", header: "immintrin.h".}
proc mm_insert_epi64*(a: m128i; i: int64; imm8: cint): m128i {.
    importc: "_mm_insert_epi64", header: "immintrin.h".}
proc mm_insert_epi8*(a: m128i; i: cint; imm8: cint): m128i {.
    importc: "_mm_insert_epi8", header: "immintrin.h".}
proc mm_insert_pi16*(a: m64; i: cint; imm8: cint): m64 {.
    importc: "_mm_insert_pi16", header: "immintrin.h".}
proc mm_insert_ps*(a: m128; b: m128; imm8: cint): m128 {.
    importc: "_mm_insert_ps", header: "immintrin.h".}
proc mm_lddqu_si128*(mem_addr: ptr m128i): m128i {.importc: "_mm_lddqu_si128", 
    header: "immintrin.h".}
proc mm_lfence*() {.importc: "_mm_lfence", header: "immintrin.h".}
proc mm_load_pd*(mem_addr: ptr cdouble): m128d {.importc: "_mm_load_pd", 
    header: "immintrin.h".}
proc mm_load_pd1*(mem_addr: ptr cdouble): m128d {.importc: "_mm_load_pd1", 
    header: "immintrin.h".}
proc mm_load_ps*(mem_addr: ptr cfloat): m128 {.importc: "_mm_load_ps", 
    header: "immintrin.h".}
proc mm_load_ps1*(mem_addr: ptr cfloat): m128 {.importc: "_mm_load_ps1", 
    header: "immintrin.h".}
proc mm_load_sd*(mem_addr: ptr cdouble): m128d {.importc: "_mm_load_sd", 
    header: "immintrin.h".}
proc mm_load_si128*(mem_addr: ptr m128i): m128i {.importc: "_mm_load_si128", 
    header: "immintrin.h".}
proc mm_load_ss*(mem_addr: ptr cfloat): m128 {.importc: "_mm_load_ss", 
    header: "immintrin.h".}
proc mm_load1_pd*(mem_addr: ptr cdouble): m128d {.importc: "_mm_load1_pd", 
    header: "immintrin.h".}
proc mm_load1_ps*(mem_addr: ptr cfloat): m128 {.importc: "_mm_load1_ps", 
    header: "immintrin.h".}
proc mm_loaddup_pd*(mem_addr: ptr cdouble): m128d {.importc: "_mm_loaddup_pd", 
    header: "immintrin.h".}
proc mm_loadh_pd*(a: m128d; mem_addr: ptr cdouble): m128d {.
    importc: "_mm_loadh_pd", header: "immintrin.h".}
proc mm_loadh_pi*(a: m128; mem_addr: ptr m64): m128 {.importc: "_mm_loadh_pi", 
    header: "immintrin.h".}
proc mm_loadl_epi64*(mem_addr: ptr m128i): m128i {.importc: "_mm_loadl_epi64", 
    header: "immintrin.h".}
proc mm_loadl_pd*(a: m128d; mem_addr: ptr cdouble): m128d {.
    importc: "_mm_loadl_pd", header: "immintrin.h".}
proc mm_loadl_pi*(a: m128; mem_addr: ptr m64): m128 {.importc: "_mm_loadl_pi", 
    header: "immintrin.h".}
proc mm_loadr_pd*(mem_addr: ptr cdouble): m128d {.importc: "_mm_loadr_pd", 
    header: "immintrin.h".}
proc mm_loadr_ps*(mem_addr: ptr cfloat): m128 {.importc: "_mm_loadr_ps", 
    header: "immintrin.h".}
proc mm_loadu_pd*(mem_addr: ptr cdouble): m128d {.importc: "_mm_loadu_pd", 
    header: "immintrin.h".}
proc mm_loadu_ps*(mem_addr: ptr cfloat): m128 {.importc: "_mm_loadu_ps", 
    header: "immintrin.h".}
proc mm_loadu_si128*(mem_addr: ptr m128i): m128i {.importc: "_mm_loadu_si128", 
    header: "immintrin.h".}
proc mm_madd_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_madd_epi16", 
    header: "immintrin.h".}
proc mm_madd_pi16*(a: m64; b: m64): m64 {.importc: "_mm_madd_pi16", 
    header: "immintrin.h".}
proc mm_maddubs_epi16*(a: m128i; b: m128i): m128i {.
    importc: "_mm_maddubs_epi16", header: "immintrin.h".}
proc mm_maddubs_pi16*(a: m64; b: m64): m64 {.importc: "_mm_maddubs_pi16", 
    header: "immintrin.h".}
proc mm_maskmove_si64*(a: m64; mask: m64; mem_addr: cstring) {.
    importc: "_mm_maskmove_si64", header: "immintrin.h".}
proc mm_maskmoveu_si128*(a: m128i; mask: m128i; mem_addr: cstring) {.
    importc: "_mm_maskmoveu_si128", header: "immintrin.h".}
proc m_maskmovq*(a: m64; mask: m64; mem_addr: cstring) {.importc: "_m_maskmovq", 
    header: "immintrin.h".}
proc mm_max_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_max_epi16", 
    header: "immintrin.h".}
proc mm_max_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_max_epi32", 
    header: "immintrin.h".}
proc mm_max_epi8*(a: m128i; b: m128i): m128i {.importc: "_mm_max_epi8", 
    header: "immintrin.h".}
proc mm_max_epu16*(a: m128i; b: m128i): m128i {.importc: "_mm_max_epu16", 
    header: "immintrin.h".}
proc mm_max_epu32*(a: m128i; b: m128i): m128i {.importc: "_mm_max_epu32", 
    header: "immintrin.h".}
proc mm_max_epu8*(a: m128i; b: m128i): m128i {.importc: "_mm_max_epu8", 
    header: "immintrin.h".}
proc mm_max_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_max_pd", 
    header: "immintrin.h".}
proc mm_max_pi16*(a: m64; b: m64): m64 {.importc: "_mm_max_pi16", 
    header: "immintrin.h".}
proc mm_max_ps*(a: m128; b: m128): m128 {.importc: "_mm_max_ps", 
    header: "immintrin.h".}
proc mm_max_pu8*(a: m64; b: m64): m64 {.importc: "_mm_max_pu8", 
                                        header: "immintrin.h".}
proc mm_max_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_max_sd", 
    header: "immintrin.h".}
proc mm_max_ss*(a: m128; b: m128): m128 {.importc: "_mm_max_ss", 
    header: "immintrin.h".}
proc mm_mfence*() {.importc: "_mm_mfence", header: "immintrin.h".}
proc mm_min_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_min_epi16", 
    header: "immintrin.h".}
proc mm_min_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_min_epi32", 
    header: "immintrin.h".}
proc mm_min_epi8*(a: m128i; b: m128i): m128i {.importc: "_mm_min_epi8", 
    header: "immintrin.h".}
proc mm_min_epu16*(a: m128i; b: m128i): m128i {.importc: "_mm_min_epu16", 
    header: "immintrin.h".}
proc mm_min_epu32*(a: m128i; b: m128i): m128i {.importc: "_mm_min_epu32", 
    header: "immintrin.h".}
proc mm_min_epu8*(a: m128i; b: m128i): m128i {.importc: "_mm_min_epu8", 
    header: "immintrin.h".}
proc mm_min_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_min_pd", 
    header: "immintrin.h".}
proc mm_min_pi16*(a: m64; b: m64): m64 {.importc: "_mm_min_pi16", 
    header: "immintrin.h".}
proc mm_min_ps*(a: m128; b: m128): m128 {.importc: "_mm_min_ps", 
    header: "immintrin.h".}
proc mm_min_pu8*(a: m64; b: m64): m64 {.importc: "_mm_min_pu8", 
                                        header: "immintrin.h".}
proc mm_min_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_min_sd", 
    header: "immintrin.h".}
proc mm_min_ss*(a: m128; b: m128): m128 {.importc: "_mm_min_ss", 
    header: "immintrin.h".}
proc mm_minpos_epu16*(a: m128i): m128i {.importc: "_mm_minpos_epu16", 
    header: "immintrin.h".}
proc mm_move_epi64*(a: m128i): m128i {.importc: "_mm_move_epi64", 
                                       header: "immintrin.h".}
proc mm_move_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_move_sd", 
    header: "immintrin.h".}
proc mm_move_ss*(a: m128; b: m128): m128 {.importc: "_mm_move_ss", 
    header: "immintrin.h".}
proc mm_movedup_pd*(a: m128d): m128d {.importc: "_mm_movedup_pd", 
                                       header: "immintrin.h".}
proc mm_movehdup_ps*(a: m128): m128 {.importc: "_mm_movehdup_ps", 
                                      header: "immintrin.h".}
proc mm_movehl_ps*(a: m128; b: m128): m128 {.importc: "_mm_movehl_ps", 
    header: "immintrin.h".}
proc mm_moveldup_ps*(a: m128): m128 {.importc: "_mm_moveldup_ps", 
                                      header: "immintrin.h".}
proc mm_movelh_ps*(a: m128; b: m128): m128 {.importc: "_mm_movelh_ps", 
    header: "immintrin.h".}
proc mm_movemask_epi8*(a: m128i): cint {.importc: "_mm_movemask_epi8", 
    header: "immintrin.h".}
proc mm_movemask_pd*(a: m128d): cint {.importc: "_mm_movemask_pd", 
                                       header: "immintrin.h".}
proc mm_movemask_pi8*(a: m64): cint {.importc: "_mm_movemask_pi8", 
                                      header: "immintrin.h".}
proc mm_movemask_ps*(a: m128): cint {.importc: "_mm_movemask_ps", 
                                      header: "immintrin.h".}
proc mm_movepi64_pi64*(a: m128i): m64 {.importc: "_mm_movepi64_pi64", 
                                        header: "immintrin.h".}
proc mm_movpi64_epi64*(a: m64): m128i {.importc: "_mm_movpi64_epi64", 
                                        header: "immintrin.h".}
proc mm_mpsadbw_epu8*(a: m128i; b: m128i; imm8: cint): m128i {.
    importc: "_mm_mpsadbw_epu8", header: "immintrin.h".}
proc mm_mul_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_mul_epi32", 
    header: "immintrin.h".}
proc mm_mul_epu32*(a: m128i; b: m128i): m128i {.importc: "_mm_mul_epu32", 
    header: "immintrin.h".}
proc mm_mul_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_mul_pd", 
    header: "immintrin.h".}
proc mm_mul_ps*(a: m128; b: m128): m128 {.importc: "_mm_mul_ps", 
    header: "immintrin.h".}
proc mm_mul_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_mul_sd", 
    header: "immintrin.h".}
proc mm_mul_ss*(a: m128; b: m128): m128 {.importc: "_mm_mul_ss", 
    header: "immintrin.h".}
proc mm_mul_su32*(a: m64; b: m64): m64 {.importc: "_mm_mul_su32", 
    header: "immintrin.h".}
proc mm_mulhi_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_mulhi_epi16", 
    header: "immintrin.h".}
proc mm_mulhi_epu16*(a: m128i; b: m128i): m128i {.importc: "_mm_mulhi_epu16", 
    header: "immintrin.h".}
proc mm_mulhi_pi16*(a: m64; b: m64): m64 {.importc: "_mm_mulhi_pi16", 
    header: "immintrin.h".}
proc mm_mulhi_pu16*(a: m64; b: m64): m64 {.importc: "_mm_mulhi_pu16", 
    header: "immintrin.h".}
proc mm_mulhrs_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_mulhrs_epi16", 
    header: "immintrin.h".}
proc mm_mulhrs_pi16*(a: m64; b: m64): m64 {.importc: "_mm_mulhrs_pi16", 
    header: "immintrin.h".}
proc mm_mullo_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_mullo_epi16", 
    header: "immintrin.h".}
proc mm_mullo_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_mullo_epi32", 
    header: "immintrin.h".}
proc mm_mullo_pi16*(a: m64; b: m64): m64 {.importc: "_mm_mullo_pi16", 
    header: "immintrin.h".}
proc mm_or_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_or_pd", 
    header: "immintrin.h".}
proc mm_or_ps*(a: m128; b: m128): m128 {.importc: "_mm_or_ps", 
    header: "immintrin.h".}
proc mm_or_si128*(a: m128i; b: m128i): m128i {.importc: "_mm_or_si128", 
    header: "immintrin.h".}
proc mm_or_si64*(a: m64; b: m64): m64 {.importc: "_mm_or_si64", 
                                        header: "immintrin.h".}
proc mm_packs_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_packs_epi16", 
    header: "immintrin.h".}
proc mm_packs_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_packs_epi32", 
    header: "immintrin.h".}
proc mm_packs_pi16*(a: m64; b: m64): m64 {.importc: "_mm_packs_pi16", 
    header: "immintrin.h".}
proc mm_packs_pi32*(a: m64; b: m64): m64 {.importc: "_mm_packs_pi32", 
    header: "immintrin.h".}
proc mm_packs_pu16*(a: m64; b: m64): m64 {.importc: "_mm_packs_pu16", 
    header: "immintrin.h".}
proc m_packssdw*(a: m64; b: m64): m64 {.importc: "_m_packssdw", 
                                        header: "immintrin.h".}
proc m_packsswb*(a: m64; b: m64): m64 {.importc: "_m_packsswb", 
                                        header: "immintrin.h".}
proc mm_packus_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_packus_epi16", 
    header: "immintrin.h".}
proc mm_packus_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_packus_epi32", 
    header: "immintrin.h".}
proc m_packuswb*(a: m64; b: m64): m64 {.importc: "_m_packuswb", 
                                        header: "immintrin.h".}
proc m_paddb*(a: m64; b: m64): m64 {.importc: "_m_paddb", header: "immintrin.h".}
proc m_paddd*(a: m64; b: m64): m64 {.importc: "_m_paddd", header: "immintrin.h".}
proc m_paddsb*(a: m64; b: m64): m64 {.importc: "_m_paddsb", 
                                      header: "immintrin.h".}
proc m_paddsw*(a: m64; b: m64): m64 {.importc: "_m_paddsw", 
                                      header: "immintrin.h".}
proc m_paddusb*(a: m64; b: m64): m64 {.importc: "_m_paddusb", 
                                       header: "immintrin.h".}
proc m_paddusw*(a: m64; b: m64): m64 {.importc: "_m_paddusw", 
                                       header: "immintrin.h".}
proc m_paddw*(a: m64; b: m64): m64 {.importc: "_m_paddw", header: "immintrin.h".}
proc m_pand*(a: m64; b: m64): m64 {.importc: "_m_pand", header: "immintrin.h".}
proc m_pandn*(a: m64; b: m64): m64 {.importc: "_m_pandn", header: "immintrin.h".}
proc mm_pause*() {.importc: "_mm_pause", header: "immintrin.h".}
proc m_pavgb*(a: m64; b: m64): m64 {.importc: "_m_pavgb", header: "immintrin.h".}
proc m_pavgw*(a: m64; b: m64): m64 {.importc: "_m_pavgw", header: "immintrin.h".}
proc m_pcmpeqb*(a: m64; b: m64): m64 {.importc: "_m_pcmpeqb", 
                                       header: "immintrin.h".}
proc m_pcmpeqd*(a: m64; b: m64): m64 {.importc: "_m_pcmpeqd", 
                                       header: "immintrin.h".}
proc m_pcmpeqw*(a: m64; b: m64): m64 {.importc: "_m_pcmpeqw", 
                                       header: "immintrin.h".}
proc m_pcmpgtb*(a: m64; b: m64): m64 {.importc: "_m_pcmpgtb", 
                                       header: "immintrin.h".}
proc m_pcmpgtd*(a: m64; b: m64): m64 {.importc: "_m_pcmpgtd", 
                                       header: "immintrin.h".}
proc m_pcmpgtw*(a: m64; b: m64): m64 {.importc: "_m_pcmpgtw", 
                                       header: "immintrin.h".}
proc m_pextrw*(a: m64; imm8: cint): cint {.importc: "_m_pextrw", 
    header: "immintrin.h".}
proc m_pinsrw*(a: m64; i: cint; imm8: cint): m64 {.importc: "_m_pinsrw", 
    header: "immintrin.h".}
proc m_pmaddwd*(a: m64; b: m64): m64 {.importc: "_m_pmaddwd", 
                                       header: "immintrin.h".}
proc m_pmaxsw*(a: m64; b: m64): m64 {.importc: "_m_pmaxsw", 
                                      header: "immintrin.h".}
proc m_pmaxub*(a: m64; b: m64): m64 {.importc: "_m_pmaxub", 
                                      header: "immintrin.h".}
proc m_pminsw*(a: m64; b: m64): m64 {.importc: "_m_pminsw", 
                                      header: "immintrin.h".}
proc m_pminub*(a: m64; b: m64): m64 {.importc: "_m_pminub", 
                                      header: "immintrin.h".}
proc m_pmovmskb*(a: m64): cint {.importc: "_m_pmovmskb", header: "immintrin.h".}
proc m_pmulhuw*(a: m64; b: m64): m64 {.importc: "_m_pmulhuw", 
                                       header: "immintrin.h".}
proc m_pmulhw*(a: m64; b: m64): m64 {.importc: "_m_pmulhw", 
                                      header: "immintrin.h".}
proc m_pmullw*(a: m64; b: m64): m64 {.importc: "_m_pmullw", 
                                      header: "immintrin.h".}
proc m_por*(a: m64; b: m64): m64 {.importc: "_m_por", header: "immintrin.h".}
proc mm_prefetch*(p: cstring; i: cint) {.importc: "_mm_prefetch", 
    header: "immintrin.h".}
proc m_psadbw*(a: m64; b: m64): m64 {.importc: "_m_psadbw", 
                                      header: "immintrin.h".}
proc m_pshufw*(a: m64; imm8: cint): m64 {.importc: "_m_pshufw", 
    header: "immintrin.h".}
proc m_pslld*(a: m64; count: m64): m64 {.importc: "_m_pslld", 
    header: "immintrin.h".}
proc m_pslldi*(a: m64; imm8: cint): m64 {.importc: "_m_pslldi", 
    header: "immintrin.h".}
proc m_psllq*(a: m64; count: m64): m64 {.importc: "_m_psllq", 
    header: "immintrin.h".}
proc m_psllqi*(a: m64; imm8: cint): m64 {.importc: "_m_psllqi", 
    header: "immintrin.h".}
proc m_psllw*(a: m64; count: m64): m64 {.importc: "_m_psllw", 
    header: "immintrin.h".}
proc m_psllwi*(a: m64; imm8: cint): m64 {.importc: "_m_psllwi", 
    header: "immintrin.h".}
proc m_psrad*(a: m64; count: m64): m64 {.importc: "_m_psrad", 
    header: "immintrin.h".}
proc m_psradi*(a: m64; imm8: cint): m64 {.importc: "_m_psradi", 
    header: "immintrin.h".}
proc m_psraw*(a: m64; count: m64): m64 {.importc: "_m_psraw", 
    header: "immintrin.h".}
proc m_psrawi*(a: m64; imm8: cint): m64 {.importc: "_m_psrawi", 
    header: "immintrin.h".}
proc m_psrld*(a: m64; count: m64): m64 {.importc: "_m_psrld", 
    header: "immintrin.h".}
proc m_psrldi*(a: m64; imm8: cint): m64 {.importc: "_m_psrldi", 
    header: "immintrin.h".}
proc m_psrlq*(a: m64; count: m64): m64 {.importc: "_m_psrlq", 
    header: "immintrin.h".}
proc m_psrlqi*(a: m64; imm8: cint): m64 {.importc: "_m_psrlqi", 
    header: "immintrin.h".}
proc m_psrlw*(a: m64; count: m64): m64 {.importc: "_m_psrlw", 
    header: "immintrin.h".}
proc m_psrlwi*(a: m64; imm8: cint): m64 {.importc: "_m_psrlwi", 
    header: "immintrin.h".}
proc m_psubb*(a: m64; b: m64): m64 {.importc: "_m_psubb", header: "immintrin.h".}
proc m_psubd*(a: m64; b: m64): m64 {.importc: "_m_psubd", header: "immintrin.h".}
proc m_psubsb*(a: m64; b: m64): m64 {.importc: "_m_psubsb", 
                                      header: "immintrin.h".}
proc m_psubsw*(a: m64; b: m64): m64 {.importc: "_m_psubsw", 
                                      header: "immintrin.h".}
proc m_psubusb*(a: m64; b: m64): m64 {.importc: "_m_psubusb", 
                                       header: "immintrin.h".}
proc m_psubusw*(a: m64; b: m64): m64 {.importc: "_m_psubusw", 
                                       header: "immintrin.h".}
proc m_psubw*(a: m64; b: m64): m64 {.importc: "_m_psubw", header: "immintrin.h".}
proc m_punpckhbw*(a: m64; b: m64): m64 {.importc: "_m_punpckhbw", 
    header: "immintrin.h".}
proc m_punpckhdq*(a: m64; b: m64): m64 {.importc: "_m_punpckhdq", 
    header: "immintrin.h".}
proc m_punpckhwd*(a: m64; b: m64): m64 {.importc: "_m_punpckhwd", 
    header: "immintrin.h".}
proc m_punpcklbw*(a: m64; b: m64): m64 {.importc: "_m_punpcklbw", 
    header: "immintrin.h".}
proc m_punpckldq*(a: m64; b: m64): m64 {.importc: "_m_punpckldq", 
    header: "immintrin.h".}
proc m_punpcklwd*(a: m64; b: m64): m64 {.importc: "_m_punpcklwd", 
    header: "immintrin.h".}
proc m_pxor*(a: m64; b: m64): m64 {.importc: "_m_pxor", header: "immintrin.h".}
proc mm_rcp_ps*(a: m128): m128 {.importc: "_mm_rcp_ps", header: "immintrin.h".}
proc mm_rcp_ss*(a: m128): m128 {.importc: "_mm_rcp_ss", header: "immintrin.h".}
proc mm_round_pd*(a: m128d; rounding: cint): m128d {.importc: "_mm_round_pd", 
    header: "immintrin.h".}
proc mm_round_ps*(a: m128; rounding: cint): m128 {.importc: "_mm_round_ps", 
    header: "immintrin.h".}
proc mm_round_sd*(a: m128d; b: m128d; rounding: cint): m128d {.
    importc: "_mm_round_sd", header: "immintrin.h".}
proc mm_round_ss*(a: m128; b: m128; rounding: cint): m128 {.
    importc: "_mm_round_ss", header: "immintrin.h".}
proc mm_rsqrt_ps*(a: m128): m128 {.importc: "_mm_rsqrt_ps", 
                                   header: "immintrin.h".}
proc mm_rsqrt_ss*(a: m128): m128 {.importc: "_mm_rsqrt_ss", 
                                   header: "immintrin.h".}
proc mm_sad_epu8*(a: m128i; b: m128i): m128i {.importc: "_mm_sad_epu8", 
    header: "immintrin.h".}
proc mm_sad_pu8*(a: m64; b: m64): m64 {.importc: "_mm_sad_pu8", 
                                        header: "immintrin.h".}
proc mm_set_epi16*(e7: cshort; e6: cshort; e5: cshort; e4: cshort; e3: cshort; 
                   e2: cshort; e1: cshort; e0: cshort): m128i {.
    importc: "_mm_set_epi16", header: "immintrin.h".}
proc mm_set_epi32*(e3: cint; e2: cint; e1: cint; e0: cint): m128i {.
    importc: "_mm_set_epi32", header: "immintrin.h".}
proc mm_set_epi64*(e1: m64; e0: m64): m128i {.importc: "_mm_set_epi64", 
    header: "immintrin.h".}
proc mm_set_epi64x*(e1: int64; e0: int64): m128i {.importc: "_mm_set_epi64x", 
    header: "immintrin.h".}
proc mm_set_epi8*(e15: char; e14: char; e13: char; e12: char; e11: char; 
                  e10: char; e9: char; e8: char; e7: char; e6: char; e5: char; 
                  e4: char; e3: char; e2: char; e1: char; e0: char): m128i {.
    importc: "_mm_set_epi8", header: "immintrin.h".}
proc MM_SET_EXCEPTION_MASK*(a: cuint) {.importc: "_MM_SET_EXCEPTION_MASK", 
                                        header: "immintrin.h".}
proc MM_SET_EXCEPTION_STATE*(a: cuint) {.importc: "_MM_SET_EXCEPTION_STATE", 
    header: "immintrin.h".}
proc MM_SET_FLUSH_ZERO_MODE*(a: cuint) {.importc: "_MM_SET_FLUSH_ZERO_MODE", 
    header: "immintrin.h".}
proc mm_set_pd*(e1: cdouble; e0: cdouble): m128d {.importc: "_mm_set_pd", 
    header: "immintrin.h".}
proc mm_set_pd1*(a: cdouble): m128d {.importc: "_mm_set_pd1", 
                                      header: "immintrin.h".}
proc mm_set_pi16*(e3: cshort; e2: cshort; e1: cshort; e0: cshort): m64 {.
    importc: "_mm_set_pi16", header: "immintrin.h".}
proc mm_set_pi32*(e1: cint; e0: cint): m64 {.importc: "_mm_set_pi32", 
    header: "immintrin.h".}
proc mm_set_pi8*(e7: char; e6: char; e5: char; e4: char; e3: char; e2: char; 
                 e1: char; e0: char): m64 {.importc: "_mm_set_pi8", 
    header: "immintrin.h".}
proc mm_set_ps*(e3: cfloat; e2: cfloat; e1: cfloat; e0: cfloat): m128 {.
    importc: "_mm_set_ps", header: "immintrin.h".}
proc mm_set_ps1*(a: cfloat): m128 {.importc: "_mm_set_ps1", 
                                    header: "immintrin.h".}
proc MM_SET_ROUNDING_MODE*(a: cuint) {.importc: "_MM_SET_ROUNDING_MODE", 
                                       header: "immintrin.h".}
proc mm_set_sd*(a: cdouble): m128d {.importc: "_mm_set_sd", 
                                     header: "immintrin.h".}
proc mm_set_ss*(a: cfloat): m128 {.importc: "_mm_set_ss", header: "immintrin.h".}
proc mm_set1_epi16*(a: cshort): m128i {.importc: "_mm_set1_epi16", 
                                        header: "immintrin.h".}
proc mm_set1_epi32*(a: cint): m128i {.importc: "_mm_set1_epi32", 
                                      header: "immintrin.h".}
proc mm_set1_epi64*(a: m64): m128i {.importc: "_mm_set1_epi64", 
                                     header: "immintrin.h".}
proc mm_set1_epi64x*(a: int64): m128i {.importc: "_mm_set1_epi64x", 
                                        header: "immintrin.h".}
proc mm_set1_epi8*(a: char): m128i {.importc: "_mm_set1_epi8", 
                                     header: "immintrin.h".}
proc mm_set1_pd*(a: cdouble): m128d {.importc: "_mm_set1_pd", 
                                      header: "immintrin.h".}
proc mm_set1_pi16*(a: cshort): m64 {.importc: "_mm_set1_pi16", 
                                     header: "immintrin.h".}
proc mm_set1_pi32*(a: cint): m64 {.importc: "_mm_set1_pi32", 
                                   header: "immintrin.h".}
proc mm_set1_pi8*(a: char): m64 {.importc: "_mm_set1_pi8", header: "immintrin.h".}
proc mm_set1_ps*(a: cfloat): m128 {.importc: "_mm_set1_ps", 
                                    header: "immintrin.h".}
proc mm_setcsr*(a: cuint) {.importc: "_mm_setcsr", header: "immintrin.h".}
proc mm_setr_epi16*(e7: cshort; e6: cshort; e5: cshort; e4: cshort; e3: cshort; 
                    e2: cshort; e1: cshort; e0: cshort): m128i {.
    importc: "_mm_setr_epi16", header: "immintrin.h".}
proc mm_setr_epi32*(e3: cint; e2: cint; e1: cint; e0: cint): m128i {.
    importc: "_mm_setr_epi32", header: "immintrin.h".}
proc mm_setr_epi64*(e1: m64; e0: m64): m128i {.importc: "_mm_setr_epi64", 
    header: "immintrin.h".}
proc mm_setr_epi8*(e15: char; e14: char; e13: char; e12: char; e11: char; 
                   e10: char; e9: char; e8: char; e7: char; e6: char; e5: char; 
                   e4: char; e3: char; e2: char; e1: char; e0: char): m128i {.
    importc: "_mm_setr_epi8", header: "immintrin.h".}
proc mm_setr_pd*(e1: cdouble; e0: cdouble): m128d {.importc: "_mm_setr_pd", 
    header: "immintrin.h".}
proc mm_setr_pi16*(e3: cshort; e2: cshort; e1: cshort; e0: cshort): m64 {.
    importc: "_mm_setr_pi16", header: "immintrin.h".}
proc mm_setr_pi32*(e1: cint; e0: cint): m64 {.importc: "_mm_setr_pi32", 
    header: "immintrin.h".}
proc mm_setr_pi8*(e7: char; e6: char; e5: char; e4: char; e3: char; e2: char; 
                  e1: char; e0: char): m64 {.importc: "_mm_setr_pi8", 
    header: "immintrin.h".}
proc mm_setr_ps*(e3: cfloat; e2: cfloat; e1: cfloat; e0: cfloat): m128 {.
    importc: "_mm_setr_ps", header: "immintrin.h".}
proc mm_setzero_pd*(): m128d {.importc: "_mm_setzero_pd", header: "immintrin.h".}
proc mm_setzero_ps*(): m128 {.importc: "_mm_setzero_ps", header: "immintrin.h".}
proc mm_setzero_si128*(): m128i {.importc: "_mm_setzero_si128", 
                                  header: "immintrin.h".}
proc mm_setzero_si64*(): m64 {.importc: "_mm_setzero_si64", 
                               header: "immintrin.h".}
proc mm_sfence*() {.importc: "_mm_sfence", header: "immintrin.h".}
proc mm_shuffle_epi32*(a: m128i; imm8: cint): m128i {.
    importc: "_mm_shuffle_epi32", header: "immintrin.h".}
proc mm_shuffle_epi8*(a: m128i; b: m128i): m128i {.importc: "_mm_shuffle_epi8", 
    header: "immintrin.h".}
proc mm_shuffle_pd*(a: m128d; b: m128d; imm8: cint): m128d {.
    importc: "_mm_shuffle_pd", header: "immintrin.h".}
proc mm_shuffle_pi16*(a: m64; imm8: cint): m64 {.importc: "_mm_shuffle_pi16", 
    header: "immintrin.h".}
proc mm_shuffle_pi8*(a: m64; b: m64): m64 {.importc: "_mm_shuffle_pi8", 
    header: "immintrin.h".}
proc mm_shuffle_ps*(a: m128; b: m128; imm8: cuint): m128 {.
    importc: "_mm_shuffle_ps", header: "immintrin.h".}
proc mm_shufflehi_epi16*(a: m128i; imm8: cint): m128i {.
    importc: "_mm_shufflehi_epi16", header: "immintrin.h".}
proc mm_shufflelo_epi16*(a: m128i; imm8: cint): m128i {.
    importc: "_mm_shufflelo_epi16", header: "immintrin.h".}
proc mm_sign_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_sign_epi16", 
    header: "immintrin.h".}
proc mm_sign_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_sign_epi32", 
    header: "immintrin.h".}
proc mm_sign_epi8*(a: m128i; b: m128i): m128i {.importc: "_mm_sign_epi8", 
    header: "immintrin.h".}
proc mm_sign_pi16*(a: m64; b: m64): m64 {.importc: "_mm_sign_pi16", 
    header: "immintrin.h".}
proc mm_sign_pi32*(a: m64; b: m64): m64 {.importc: "_mm_sign_pi32", 
    header: "immintrin.h".}
proc mm_sign_pi8*(a: m64; b: m64): m64 {.importc: "_mm_sign_pi8", 
    header: "immintrin.h".}
proc mm_sll_epi16*(a: m128i; count: m128i): m128i {.importc: "_mm_sll_epi16", 
    header: "immintrin.h".}
proc mm_sll_epi32*(a: m128i; count: m128i): m128i {.importc: "_mm_sll_epi32", 
    header: "immintrin.h".}
proc mm_sll_epi64*(a: m128i; count: m128i): m128i {.importc: "_mm_sll_epi64", 
    header: "immintrin.h".}
proc mm_sll_pi16*(a: m64; count: m64): m64 {.importc: "_mm_sll_pi16", 
    header: "immintrin.h".}
proc mm_sll_pi32*(a: m64; count: m64): m64 {.importc: "_mm_sll_pi32", 
    header: "immintrin.h".}
proc mm_sll_si64*(a: m64; count: m64): m64 {.importc: "_mm_sll_si64", 
    header: "immintrin.h".}
proc mm_slli_epi16*(a: m128i; imm8: cint): m128i {.importc: "_mm_slli_epi16", 
    header: "immintrin.h".}
proc mm_slli_epi32*(a: m128i; imm8: cint): m128i {.importc: "_mm_slli_epi32", 
    header: "immintrin.h".}
proc mm_slli_epi64*(a: m128i; imm8: cint): m128i {.importc: "_mm_slli_epi64", 
    header: "immintrin.h".}
proc mm_slli_pi16*(a: m64; imm8: cint): m64 {.importc: "_mm_slli_pi16", 
    header: "immintrin.h".}
proc mm_slli_pi32*(a: m64; imm8: cint): m64 {.importc: "_mm_slli_pi32", 
    header: "immintrin.h".}
proc mm_slli_si128*(a: m128i; imm8: cint): m128i {.importc: "_mm_slli_si128", 
    header: "immintrin.h".}
proc mm_slli_si64*(a: m64; imm8: cint): m64 {.importc: "_mm_slli_si64", 
    header: "immintrin.h".}
proc mm_sqrt_pd*(a: m128d): m128d {.importc: "_mm_sqrt_pd", 
                                    header: "immintrin.h".}
proc mm_sqrt_ps*(a: m128): m128 {.importc: "_mm_sqrt_ps", header: "immintrin.h".}
proc mm_sqrt_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_sqrt_sd", 
    header: "immintrin.h".}
proc mm_sqrt_ss*(a: m128): m128 {.importc: "_mm_sqrt_ss", header: "immintrin.h".}
proc mm_sra_epi16*(a: m128i; count: m128i): m128i {.importc: "_mm_sra_epi16", 
    header: "immintrin.h".}
proc mm_sra_epi32*(a: m128i; count: m128i): m128i {.importc: "_mm_sra_epi32", 
    header: "immintrin.h".}
proc mm_sra_pi16*(a: m64; count: m64): m64 {.importc: "_mm_sra_pi16", 
    header: "immintrin.h".}
proc mm_sra_pi32*(a: m64; count: m64): m64 {.importc: "_mm_sra_pi32", 
    header: "immintrin.h".}
proc mm_srai_epi16*(a: m128i; imm8: cint): m128i {.importc: "_mm_srai_epi16", 
    header: "immintrin.h".}
proc mm_srai_epi32*(a: m128i; imm8: cint): m128i {.importc: "_mm_srai_epi32", 
    header: "immintrin.h".}
proc mm_srai_pi16*(a: m64; imm8: cint): m64 {.importc: "_mm_srai_pi16", 
    header: "immintrin.h".}
proc mm_srai_pi32*(a: m64; imm8: cint): m64 {.importc: "_mm_srai_pi32", 
    header: "immintrin.h".}
proc mm_srl_epi16*(a: m128i; count: m128i): m128i {.importc: "_mm_srl_epi16", 
    header: "immintrin.h".}
proc mm_srl_epi32*(a: m128i; count: m128i): m128i {.importc: "_mm_srl_epi32", 
    header: "immintrin.h".}
proc mm_srl_epi64*(a: m128i; count: m128i): m128i {.importc: "_mm_srl_epi64", 
    header: "immintrin.h".}
proc mm_srl_pi16*(a: m64; count: m64): m64 {.importc: "_mm_srl_pi16", 
    header: "immintrin.h".}
proc mm_srl_pi32*(a: m64; count: m64): m64 {.importc: "_mm_srl_pi32", 
    header: "immintrin.h".}
proc mm_srl_si64*(a: m64; count: m64): m64 {.importc: "_mm_srl_si64", 
    header: "immintrin.h".}
proc mm_srli_epi16*(a: m128i; imm8: cint): m128i {.importc: "_mm_srli_epi16", 
    header: "immintrin.h".}
proc mm_srli_epi32*(a: m128i; imm8: cint): m128i {.importc: "_mm_srli_epi32", 
    header: "immintrin.h".}
proc mm_srli_epi64*(a: m128i; imm8: cint): m128i {.importc: "_mm_srli_epi64", 
    header: "immintrin.h".}
proc mm_srli_pi16*(a: m64; imm8: cint): m64 {.importc: "_mm_srli_pi16", 
    header: "immintrin.h".}
proc mm_srli_pi32*(a: m64; imm8: cint): m64 {.importc: "_mm_srli_pi32", 
    header: "immintrin.h".}
proc mm_srli_si128*(a: m128i; imm8: cint): m128i {.importc: "_mm_srli_si128", 
    header: "immintrin.h".}
proc mm_srli_si64*(a: m64; imm8: cint): m64 {.importc: "_mm_srli_si64", 
    header: "immintrin.h".}
proc mm_store_pd*(mem_addr: ptr cdouble; a: m128d) {.importc: "_mm_store_pd", 
    header: "immintrin.h".}
proc mm_store_pd1*(mem_addr: ptr cdouble; a: m128d) {.importc: "_mm_store_pd1", 
    header: "immintrin.h".}
proc mm_store_ps*(mem_addr: ptr cfloat; a: m128) {.importc: "_mm_store_ps", 
    header: "immintrin.h".}
proc mm_store_ps1*(mem_addr: ptr cfloat; a: m128) {.importc: "_mm_store_ps1", 
    header: "immintrin.h".}
proc mm_store_sd*(mem_addr: ptr cdouble; a: m128d) {.importc: "_mm_store_sd", 
    header: "immintrin.h".}
proc mm_store_si128*(mem_addr: ptr m128i; a: m128i) {.
    importc: "_mm_store_si128", header: "immintrin.h".}
proc mm_store_ss*(mem_addr: ptr cfloat; a: m128) {.importc: "_mm_store_ss", 
    header: "immintrin.h".}
proc mm_store1_pd*(mem_addr: ptr cdouble; a: m128d) {.importc: "_mm_store1_pd", 
    header: "immintrin.h".}
proc mm_store1_ps*(mem_addr: ptr cfloat; a: m128) {.importc: "_mm_store1_ps", 
    header: "immintrin.h".}
proc mm_storeh_pd*(mem_addr: ptr cdouble; a: m128d) {.importc: "_mm_storeh_pd", 
    header: "immintrin.h".}
proc mm_storeh_pi*(mem_addr: ptr m64; a: m128) {.importc: "_mm_storeh_pi", 
    header: "immintrin.h".}
proc mm_storel_epi64*(mem_addr: ptr m128i; a: m128i) {.
    importc: "_mm_storel_epi64", header: "immintrin.h".}
proc mm_storel_pd*(mem_addr: ptr cdouble; a: m128d) {.importc: "_mm_storel_pd", 
    header: "immintrin.h".}
proc mm_storel_pi*(mem_addr: ptr m64; a: m128) {.importc: "_mm_storel_pi", 
    header: "immintrin.h".}
proc mm_storer_pd*(mem_addr: ptr cdouble; a: m128d) {.importc: "_mm_storer_pd", 
    header: "immintrin.h".}
proc mm_storer_ps*(mem_addr: ptr cfloat; a: m128) {.importc: "_mm_storer_ps", 
    header: "immintrin.h".}
proc mm_storeu_pd*(mem_addr: ptr cdouble; a: m128d) {.importc: "_mm_storeu_pd", 
    header: "immintrin.h".}
proc mm_storeu_ps*(mem_addr: ptr cfloat; a: m128) {.importc: "_mm_storeu_ps", 
    header: "immintrin.h".}
proc mm_storeu_si128*(mem_addr: ptr m128i; a: m128i) {.
    importc: "_mm_storeu_si128", header: "immintrin.h".}
proc mm_stream_load_si128*(mem_addr: ptr m128i): m128i {.
    importc: "_mm_stream_load_si128", header: "immintrin.h".}
proc mm_stream_pd*(mem_addr: ptr cdouble; a: m128d) {.importc: "_mm_stream_pd", 
    header: "immintrin.h".}
proc mm_stream_pi*(mem_addr: ptr m64; a: m64) {.importc: "_mm_stream_pi", 
    header: "immintrin.h".}
proc mm_stream_ps*(mem_addr: ptr cfloat; a: m128) {.importc: "_mm_stream_ps", 
    header: "immintrin.h".}
proc mm_stream_si128*(mem_addr: ptr m128i; a: m128i) {.
    importc: "_mm_stream_si128", header: "immintrin.h".}
proc mm_stream_si32*(mem_addr: ptr cint; a: cint) {.importc: "_mm_stream_si32", 
    header: "immintrin.h".}
proc mm_stream_si64*(mem_addr: ptr int64; a: int64) {.
    importc: "_mm_stream_si64", header: "immintrin.h".}
proc mm_sub_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_sub_epi16", 
    header: "immintrin.h".}
proc mm_sub_epi32*(a: m128i; b: m128i): m128i {.importc: "_mm_sub_epi32", 
    header: "immintrin.h".}
proc mm_sub_epi64*(a: m128i; b: m128i): m128i {.importc: "_mm_sub_epi64", 
    header: "immintrin.h".}
proc mm_sub_epi8*(a: m128i; b: m128i): m128i {.importc: "_mm_sub_epi8", 
    header: "immintrin.h".}
proc mm_sub_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_sub_pd", 
    header: "immintrin.h".}
proc mm_sub_pi16*(a: m64; b: m64): m64 {.importc: "_mm_sub_pi16", 
    header: "immintrin.h".}
proc mm_sub_pi32*(a: m64; b: m64): m64 {.importc: "_mm_sub_pi32", 
    header: "immintrin.h".}
proc mm_sub_pi8*(a: m64; b: m64): m64 {.importc: "_mm_sub_pi8", 
                                        header: "immintrin.h".}
proc mm_sub_ps*(a: m128; b: m128): m128 {.importc: "_mm_sub_ps", 
    header: "immintrin.h".}
proc mm_sub_sd*(a: m128d; b: m128d): m128d {.importc: "_mm_sub_sd", 
    header: "immintrin.h".}
proc mm_sub_si64*(a: m64; b: m64): m64 {.importc: "_mm_sub_si64", 
    header: "immintrin.h".}
proc mm_sub_ss*(a: m128; b: m128): m128 {.importc: "_mm_sub_ss", 
    header: "immintrin.h".}
proc mm_subs_epi16*(a: m128i; b: m128i): m128i {.importc: "_mm_subs_epi16", 
    header: "immintrin.h".}
proc mm_subs_epi8*(a: m128i; b: m128i): m128i {.importc: "_mm_subs_epi8", 
    header: "immintrin.h".}
proc mm_subs_epu16*(a: m128i; b: m128i): m128i {.importc: "_mm_subs_epu16", 
    header: "immintrin.h".}
proc mm_subs_epu8*(a: m128i; b: m128i): m128i {.importc: "_mm_subs_epu8", 
    header: "immintrin.h".}
proc mm_subs_pi16*(a: m64; b: m64): m64 {.importc: "_mm_subs_pi16", 
    header: "immintrin.h".}
proc mm_subs_pi8*(a: m64; b: m64): m64 {.importc: "_mm_subs_pi8", 
    header: "immintrin.h".}
proc mm_subs_pu16*(a: m64; b: m64): m64 {.importc: "_mm_subs_pu16", 
    header: "immintrin.h".}
proc mm_subs_pu8*(a: m64; b: m64): m64 {.importc: "_mm_subs_pu8", 
    header: "immintrin.h".}
proc mm_test_all_ones*(a: m128i): cint {.importc: "_mm_test_all_ones", 
    header: "immintrin.h".}
proc mm_test_all_zeros*(a: m128i; mask: m128i): cint {.
    importc: "_mm_test_all_zeros", header: "immintrin.h".}
proc mm_test_mix_ones_zeros*(a: m128i; mask: m128i): cint {.
    importc: "_mm_test_mix_ones_zeros", header: "immintrin.h".}
proc mm_testc_si128*(a: m128i; b: m128i): cint {.importc: "_mm_testc_si128", 
    header: "immintrin.h".}
proc mm_testnzc_si128*(a: m128i; b: m128i): cint {.importc: "_mm_testnzc_si128", 
    header: "immintrin.h".}
proc mm_testz_si128*(a: m128i; b: m128i): cint {.importc: "_mm_testz_si128", 
    header: "immintrin.h".}
proc m_to_int*(a: m64): cint {.importc: "_m_to_int", header: "immintrin.h".}
proc m_to_int64*(a: m64): int64 {.importc: "_m_to_int64", header: "immintrin.h".}
proc MM_TRANSPOSE4_PS*(row0: m128; row1: m128; row2: m128; row3: m128) {.
    importc: "_MM_TRANSPOSE4_PS", header: "immintrin.h".}
proc mm_ucomieq_sd*(a: m128d; b: m128d): cint {.importc: "_mm_ucomieq_sd", 
    header: "immintrin.h".}
proc mm_ucomieq_ss*(a: m128; b: m128): cint {.importc: "_mm_ucomieq_ss", 
    header: "immintrin.h".}
proc mm_ucomige_sd*(a: m128d; b: m128d): cint {.importc: "_mm_ucomige_sd", 
    header: "immintrin.h".}
proc mm_ucomige_ss*(a: m128; b: m128): cint {.importc: "_mm_ucomige_ss", 
    header: "immintrin.h".}
proc mm_ucomigt_sd*(a: m128d; b: m128d): cint {.importc: "_mm_ucomigt_sd", 
    header: "immintrin.h".}
proc mm_ucomigt_ss*(a: m128; b: m128): cint {.importc: "_mm_ucomigt_ss", 
    header: "immintrin.h".}
proc mm_ucomile_sd*(a: m128d; b: m128d): cint {.importc: "_mm_ucomile_sd", 
    header: "immintrin.h".}
proc mm_ucomile_ss*(a: m128; b: m128): cint {.importc: "_mm_ucomile_ss", 
    header: "immintrin.h".}
proc mm_ucomilt_sd*(a: m128d; b: m128d): cint {.importc: "_mm_ucomilt_sd", 
    header: "immintrin.h".}
proc mm_ucomilt_ss*(a: m128; b: m128): cint {.importc: "_mm_ucomilt_ss", 
    header: "immintrin.h".}
proc mm_ucomineq_sd*(a: m128d; b: m128d): cint {.importc: "_mm_ucomineq_sd", 
    header: "immintrin.h".}
proc mm_ucomineq_ss*(a: m128; b: m128): cint {.importc: "_mm_ucomineq_ss", 
    header: "immintrin.h".}
proc mm_unpackhi_epi16*(a: m128i; b: m128i): m128i {.
    importc: "_mm_unpackhi_epi16", header: "immintrin.h".}
proc mm_unpackhi_epi32*(a: m128i; b: m128i): m128i {.
    importc: "_mm_unpackhi_epi32", header: "immintrin.h".}
proc mm_unpackhi_epi64*(a: m128i; b: m128i): m128i {.
    importc: "_mm_unpackhi_epi64", header: "immintrin.h".}
proc mm_unpackhi_epi8*(a: m128i; b: m128i): m128i {.
    importc: "_mm_unpackhi_epi8", header: "immintrin.h".}
proc mm_unpackhi_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_unpackhi_pd", 
    header: "immintrin.h".}
proc mm_unpackhi_pi16*(a: m64; b: m64): m64 {.importc: "_mm_unpackhi_pi16", 
    header: "immintrin.h".}
proc mm_unpackhi_pi32*(a: m64; b: m64): m64 {.importc: "_mm_unpackhi_pi32", 
    header: "immintrin.h".}
proc mm_unpackhi_pi8*(a: m64; b: m64): m64 {.importc: "_mm_unpackhi_pi8", 
    header: "immintrin.h".}
proc mm_unpackhi_ps*(a: m128; b: m128): m128 {.importc: "_mm_unpackhi_ps", 
    header: "immintrin.h".}
proc mm_unpacklo_epi16*(a: m128i; b: m128i): m128i {.
    importc: "_mm_unpacklo_epi16", header: "immintrin.h".}
proc mm_unpacklo_epi32*(a: m128i; b: m128i): m128i {.
    importc: "_mm_unpacklo_epi32", header: "immintrin.h".}
proc mm_unpacklo_epi64*(a: m128i; b: m128i): m128i {.
    importc: "_mm_unpacklo_epi64", header: "immintrin.h".}
proc mm_unpacklo_epi8*(a: m128i; b: m128i): m128i {.
    importc: "_mm_unpacklo_epi8", header: "immintrin.h".}
proc mm_unpacklo_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_unpacklo_pd", 
    header: "immintrin.h".}
proc mm_unpacklo_pi16*(a: m64; b: m64): m64 {.importc: "_mm_unpacklo_pi16", 
    header: "immintrin.h".}
proc mm_unpacklo_pi32*(a: m64; b: m64): m64 {.importc: "_mm_unpacklo_pi32", 
    header: "immintrin.h".}
proc mm_unpacklo_pi8*(a: m64; b: m64): m64 {.importc: "_mm_unpacklo_pi8", 
    header: "immintrin.h".}
proc mm_unpacklo_ps*(a: m128; b: m128): m128 {.importc: "_mm_unpacklo_ps", 
    header: "immintrin.h".}
proc mm_xor_pd*(a: m128d; b: m128d): m128d {.importc: "_mm_xor_pd", 
    header: "immintrin.h".}
proc mm_xor_ps*(a: m128; b: m128): m128 {.importc: "_mm_xor_ps", 
    header: "immintrin.h".}
proc mm_xor_si128*(a: m128i; b: m128i): m128i {.importc: "_mm_xor_si128", 
    header: "immintrin.h".}
proc mm_xor_si64*(a: m64; b: m64): m64 {.importc: "_mm_xor_si64", 
    header: "immintrin.h".}