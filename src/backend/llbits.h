#ifndef _QEX_BACKEND_LLBITS_H_
#define _QEX_BACKEND_LLBITS_H_
#include<stdint.h>
// Prepare to break the strict aliasing rule.
typedef uint32_t __attribute__((__may_alias__,__aligned__(4))) RegisterWord;
typedef struct MemoryWord1 {RegisterWord a[1];} __attribute__((__may_alias__,__aligned__(4))) MemoryWord1;
typedef struct MemoryWord2 {RegisterWord a[2];} __attribute__((__may_alias__,__aligned__(8))) MemoryWord2;
typedef struct MemoryWord4 {RegisterWord a[4];} __attribute__((__may_alias__,__aligned__(16))) MemoryWord4;
typedef struct MemoryWord8 {RegisterWord a[8];} __attribute__((__may_alias__,__aligned__(32))) MemoryWord8;
typedef struct MemoryWord16 {RegisterWord a[16];} __attribute__((__may_alias__,__aligned__(64))) MemoryWord16;
typedef struct MemoryWord32 {RegisterWord a[32];} __attribute__((__may_alias__,__aligned__(128))) MemoryWord32;
typedef struct ShortVectorFloat1   {float   a[1] ;} __attribute__((__aligned__(1*4)))  ShortVectorFloat1  ;
typedef struct ShortVectorFloat2   {float   a[2] ;} __attribute__((__aligned__(2*4)))  ShortVectorFloat2  ;
typedef struct ShortVectorFloat4   {float   a[4] ;} __attribute__((__aligned__(4*4)))  ShortVectorFloat4  ;
typedef struct ShortVectorFloat8   {float   a[8] ;} __attribute__((__aligned__(8*4)))  ShortVectorFloat8  ;
typedef struct ShortVectorFloat16  {float   a[16];} __attribute__((__aligned__(16*4))) ShortVectorFloat16 ;
typedef struct ShortVectorFloat32  {float   a[32];} __attribute__((__aligned__(32*4))) ShortVectorFloat32 ;
typedef struct ShortVectorFloat64  {float   a[64];} __attribute__((__aligned__(64*4))) ShortVectorFloat64 ;
typedef struct ShortVectorDouble1  {double  a[1] ;} __attribute__((__aligned__(1*8)))  ShortVectorDouble1 ;
typedef struct ShortVectorDouble2  {double  a[2] ;} __attribute__((__aligned__(2*8)))  ShortVectorDouble2 ;
typedef struct ShortVectorDouble4  {double  a[4] ;} __attribute__((__aligned__(4*8)))  ShortVectorDouble4 ;
typedef struct ShortVectorDouble8  {double  a[8] ;} __attribute__((__aligned__(8*8)))  ShortVectorDouble8 ;
typedef struct ShortVectorDouble16 {double  a[16];} __attribute__((__aligned__(16*8))) ShortVectorDouble16;
typedef struct ShortVectorDouble32 {double  a[32];} __attribute__((__aligned__(32*8))) ShortVectorDouble32;
typedef struct ShortVectorInt31t1  {int32_t a[1] ;} __attribute__((__aligned__(1*4)))  ShortVectorInt31t1 ;
typedef struct ShortVectorInt31t2  {int32_t a[2] ;} __attribute__((__aligned__(2*4)))  ShortVectorInt31t2 ;
typedef struct ShortVectorInt31t4  {int32_t a[4] ;} __attribute__((__aligned__(4*4)))  ShortVectorInt31t4 ;
typedef struct ShortVectorInt31t8  {int32_t a[8] ;} __attribute__((__aligned__(8*4)))  ShortVectorInt31t8 ;
typedef struct ShortVectorInt31t16 {int32_t a[16];} __attribute__((__aligned__(16*4))) ShortVectorInt31t16;
typedef struct ShortVectorInt31t32 {int32_t a[32];} __attribute__((__aligned__(32*4))) ShortVectorInt31t32;
typedef struct ShortVectorInt31t64 {int32_t a[64];} __attribute__((__aligned__(64*4))) ShortVectorInt31t64;
typedef struct ShortVectorInt63t1  {int64_t a[1] ;} __attribute__((__aligned__(1*8)))  ShortVectorInt63t1 ;
typedef struct ShortVectorInt63t2  {int64_t a[2] ;} __attribute__((__aligned__(2*8)))  ShortVectorInt63t2 ;
typedef struct ShortVectorInt63t4  {int64_t a[4] ;} __attribute__((__aligned__(4*8)))  ShortVectorInt63t4 ;
typedef struct ShortVectorInt63t8  {int64_t a[8] ;} __attribute__((__aligned__(8*8)))  ShortVectorInt63t8 ;
typedef struct ShortVectorInt63t16 {int64_t a[16];} __attribute__((__aligned__(16*8))) ShortVectorInt63t16;
typedef struct ShortVectorInt63t32 {int64_t a[32];} __attribute__((__aligned__(32*8))) ShortVectorInt63t32;
#endif//_QEX_BACKEND_LLBITS_H_
