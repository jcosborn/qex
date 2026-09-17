## Checkpoint and latent-angle adapters for the NNFT application.
import qex
import io/[arrays, arrayfields]
import field/matrixFields
import ../../nn
import flow
import std/json

proc loadNnft*[T: SomeFloat](dir: string; precision: typedesc[T]): NnftParams[T] =
  ## Leaves param_0, param_1, ... in stage order: for each convolution its bias
  ## [out] then weights [out,in,k0,k1], and a scale [nnftCoefs,1,1] closing the
  ## stage. The shapes fix the widths, kernels and stage count.
  let man = readManifest(dir)
  if man["version"].getInt != 1:
    raise newException(ValueError, "unsupported NNFT checkpoint version")
  let leaves = man["arrays"]
  proc spec(i: int): JsonNode =
    if not leaves.hasKey("param_" & $i):
      raise newException(ValueError, "NNFT checkpoint ends inside a stage at param_" & $i)
    result = leaves["param_" & $i]
    if result["dtype"].getStr != "float32":
      raise newException(ValueError, "NNFT checkpoint arrays require float32 storage")
  var i = 0
  while leaves.hasKey("param_" & $i):
    var net: NnftNet[T]
    while true:
      let shape = arrayShape(spec(i))
      if shape.len == 3:
        if shape != @[nnftCoefs,1,1]:
          raise newException(ValueError, "NNFT scale param_" & $i & " must have shape [" & $nnftCoefs & ",1,1]")
        net.scale = arrays.readArray[T](dir,spec(i))
        inc i
        break
      let w = arrayShape(spec(i+1))
      if shape.len != 1 or w.len != 4 or w[0] != shape[0]:
        raise newException(ValueError, "NNFT param_" & $i & " must be a bias [out] followed by weights [out,in,k0,k1]")
      net.layers.add convParams(w[1],w[0],w[2..3],arrays.readArray[T](dir,spec(i+1)),arrays.readArray[T](dir,spec(i)))
      i += 2
    net.requireNet
    result.add net
  if result.len == 0:
    raise newException(ValueError, "NNFT checkpoint has no parameters")

proc loadNnft*(dir: string): NnftParams[float32] = loadNnft(dir,float32)

proc loadAngles*(g: seq[DLatticeColorMatrixV]; dir: string; spec: JsonNode) =
  ## Global C-order [direction,row,column] angles theta -> V_d = exp(i theta_d).
  g.requireNnftGauge
  let lo = g[0].l
  let theta = @[lo.RealD(),lo.RealD()]
  loadFields(theta,dir,spec)
  threads:
    for d in 0..1: expi(g[d],theta[d])

proc loadLatent*(g: seq[DLatticeColorMatrixV]; dir,name: string) =
  let man = readManifest(dir)
  if not man["arrays"].hasKey(name):
    raise newException(ValueError, "latent array is absent from checkpoint manifest: " & name)
  loadAngles(g,dir,man["arrays"][name])
