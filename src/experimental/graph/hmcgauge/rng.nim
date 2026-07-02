import qex
import strutils

type
  RngKind* = enum
    rkPhilox4x64, rkThreefry4x64, rkMrg32k3a

template withRng*(kind: RngKind; R: untyped; body: untyped) =
  case kind
  of rkPhilox4x64:
    block:
      type R {.inject.} = Philox4x64
      body
  of rkThreefry4x64:
    block:
      type R {.inject.} = Threefry4x64
      body
  of rkMrg32k3a:
    block:
      type R {.inject.} = MRG32k3a
      body

func `$`*(kind: RngKind): string =
  case kind
  of rkPhilox4x64: "Philox4x64"
  of rkThreefry4x64: "Threefry4x64"
  of rkMrg32k3a: "MRG32k3a"

func parseRngKind*(name: string): RngKind =
  case name.toLowerAscii
  of "philox4x64": rkPhilox4x64
  of "threefry4x64": rkThreefry4x64
  of "mrg32k3a": rkMrg32k3a
  else:
    raise newException(ValueError, "unknown RNG: " & name & " (expected Philox4x64, Threefry4x64, or MRG32k3a)")
