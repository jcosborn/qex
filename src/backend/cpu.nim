import base/threading

#type
#  Buffer*[T] = object
#    data: ptr UncheckedArray[T]
#    n: int

#proc init[T](b: Buffer[T], n: int) =
#  let p = allocShared(n*sizeof(T))
#  b.data = cast[type b.data](p)
#  b.n = n

#proc free[T](b: Buffer[T]) =
#  free b.data

template onGpu*(x:untyped) = threads: x
template onGpu*(n,x:untyped) = threads: x
template onGpu*(n,t,x:untyped) = threads: x
