import tensorwrap
export tensorwrap

type
  ColorType* = object
  Color*[T] = TensorObj[ColorType,T]
  Color2*[T] = Color[T]
  Color3*[T] = Color[T]
  Color4*[T] = Color[T]

#[
template getNc*[K,T](x: SomeTensor[K,T]): int =
  when K is ColorType:
    when T is Mat1:
      x[].nrows
    elif T is Vec1:
      x[].len
    else:
      static:
        echo "error: unknown Nc"
        echo x.repr
        echo type(x).name
        qexExit 1
      0
  else:
    getNc(x[])
]#

template getNc*[K,T](x: typedesc[SomeTensor[K,T]]): int =
  getNc(x[])
template getNc*[K,T](x: SomeTensor[K,T]): int =
  getNc(x[])

template getNc*[T](x: typedesc[SomeTensor[ColorType,T]]): int =
  when T is Mat1:
    x[].nrows
  elif T is Vec1:
    x[].len
  else:
    static:
      echo "error: unknown Nc"
      echo x.repr
      echo type(x).name
      qexExit 1
    0
template getNc*[T](x: SomeTensor[ColorType,T]): int =
  when T is Mat1:
    x[].nrows
  elif T is Vec1:
    x[].len
  else:
    static:
      echo "error: unknown Nc"
      echo x.repr
      echo type(x).name
      qexExit 1
    0
