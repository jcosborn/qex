## FNV-1a fingerprints.  One implementation for the frozen-rational hash
## (ops/zolotarev) and the checkpoint payload hash (hmc/trajectory); both are
## on-disk contracts, so the byte order below is fixed: little end first.

const
  fnvBasis* = 0xcbf29ce484222325'u64
  fnvPrime = 0x100000001b3'u64

func fnv1a*(h: uint64, v: uint64): uint64 =
  ## Fold the 8 bytes of v into h, least significant byte first.
  result = h
  var x = v
  for _ in 0..7:
    result = (result xor (x and 0xff'u64))*fnvPrime
    x = x shr 8

func fnv1a*(s: string): uint64 =
  ## Hash of a byte string from the basis offset.
  result = fnvBasis
  for ch in s:
    result = (result xor uint64(ord(ch)))*fnvPrime
