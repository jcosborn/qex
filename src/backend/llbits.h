#ifndef _CUDANIM_LLBITS_H_
#define _CUDANIM_LLBITS_H_
#include<stdint.h>
// Prepare to break the strict aliasing rule.
typedef uint32_t __attribute__((__may_alias__,__aligned__(4))) RegisterWord;
typedef struct MemoryWord1 {RegisterWord a[1];} __attribute__((__may_alias__,__aligned__(4))) MemoryWord1;
typedef struct MemoryWord2 {RegisterWord a[2];} __attribute__((__may_alias__,__aligned__(8))) MemoryWord2;
typedef struct MemoryWord4 {RegisterWord a[4];} __attribute__((__may_alias__,__aligned__(16))) MemoryWord4;
typedef struct MemoryWord8 {RegisterWord a[8];} __attribute__((__may_alias__,__aligned__(32))) MemoryWord8;
typedef struct MemoryWord16 {RegisterWord a[16];} __attribute__((__may_alias__,__aligned__(64))) MemoryWord16;
typedef struct MemoryWord32 {RegisterWord a[32];} __attribute__((__may_alias__,__aligned__(128))) MemoryWord32;
#endif//_CUDANIM_LLBITS_H_
