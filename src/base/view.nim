type
  View*[T] = object
    view: ptr T

template view*[T](x: typedesc[T]): typedesc = View[T]
template view*[T](x: var T): auto = View[T](view: addr x)

template `[]`*[T](x: typedesc[View[T]]): typedesc = T
template `[]`*[T](x: View[T]): T = x.view[]

template toView*(x: var SomeNumber): auto = view(x)
template toView*(x: var string): auto = view(x)
