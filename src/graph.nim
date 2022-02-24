## Runtime symbolic computational graph.

type
  Graph = ref object of RootObj
    tag: GraphTags
    val: GraphValueRef
    grad: GraphValueRef
    initGrad: proc(g:Graph)
    str: string
    case isop: bool
    of true:
      args: seq[Graph]
      run: proc(g:Graph)
      back: proc(g:Graph)
      refs: int
    of false:
      discard
  GraphTag = enum
    gtConst, gtRun, gtDF, gtGrad
  GraphTags = set[GraphTag]
  GraphValueRef = ref object of RootObj
  GraphValue[T] = ref object of GraphValueRef
    v: T
  GraphNode[T] = ref object of Graph

proc countRefs(g:Graph) =
  ## Count only go in to the children of nodes at the first encounter.
  proc clear(g:Graph) =
    if g.isop:
      g.refs = 0
      for x in g.args:
        x.clear
  proc count(g:Graph) =
    if g.isop:
      g.refs.inc
      if g.refs==1:
        for x in g.args:
          x.count
  g.clear
  g.count

proc `$`(g:Graph):string =
  g.countRefs
  var id = 0
  proc go(g:Graph):string =
    if g.isop:
      if g.refs<=0:
        result = "#" & $(-g.refs) & "#"
        return
      if g.refs>1:
        result = "#" & $id & "=("
        g.refs = -id
        inc id
      else:
        result = "("
      result &= g.str
      for x in g.args:
        result &= " " & x.go
      result &= ")"
    else:
      result = g.str
  result = g.go

proc newVar[T](x:T, s="$V"):GraphNode[T] =
  GraphNode[T](val:GraphValue[T](v:x), str:s,
    initGrad: (proc(g:Graph) = g.grad = GraphValue[T](v:1.T)))

proc newConst[T](x:T, s="$C"):GraphNode[T] =
  GraphNode[T](tag: {gtConst}, val:GraphValue[T](v:x), str:s & "|" & $x & "|",
    initGrad: (proc(g:Graph) = g.grad = GraphValue[T](v:1.T)))

proc eval(g:Graph) =
  if g.isop and gtRun notin g.tag:
    g.tag.incl gtRun
    for x in g.args:
      x.eval
    g.run(g)

proc clearGrad(g:Graph) =
  g.grad = nil
  g.tag.excl gtDF
  g.tag.excl gtGrad
  if g.isop:
    for x in g.args:
      x.clearGrad

proc evalGrad(g:Graph) =
  if g.isop and gtRun notin g.tag:
    g.eval
  if gtDF notin g.tag:
    g.clearGrad
    g.countRefs
    g.tag.incl gtDF
    g.initGrad(g)
    proc go(g:Graph) =
      g.tag.incl gtGrad
      if g.isop:
        g.refs.dec
        if g.refs>0:
          return
        # Wait until the last reference of the shared nodes.
        g.back(g)
        for x in g.args:
          x.go
    g.go

proc evalGrad[G,X](g:GraphNode[G], x:GraphNode[X]):X =
  # TODO only descend in to nodes that contains x.
  g.Graph.evalGrad
  proc go(g:Graph):Graph =
    if g == x:
      return x
    elif g.isop:
      for c in g.args:
        if c.go == x:
          return x
    g
  if g.go == x:
    GraphValue[X](x.grad).v
  else:
    X 0

proc eval[T](g:GraphNode[T]):T =
  g.Graph.eval
  GraphValue[T](g.val).v

proc `+`[X,Y](x:GraphNode[X], y:GraphNode[Y]):auto =
  type R = type(GraphValue[X](x.val).v+GraphValue[Y](y.val).v)
  GraphNode[R](isop:true, str:"+", args: @[x.Graph,y],
    run: (proc(g:Graph) =
      echo "# Run: ",g.args[0].str," + ",g.args[1].str
      let v = GraphValue[R](v:GraphValue[X](g.args[0].val).v+GraphValue[Y](g.args[1].val).v)
      g.val = v
      g.str &= "(=" & $v.v & ")"),
    initGrad: (proc(g:Graph) = g.grad = GraphValue[R](v:1.R)),
    back: (proc(g:Graph) =
      echo "# Back: ",x.str," * ",y.str
      if gtConst notin g.args[0].tag:
        let t = GraphValue[R](g.grad).v.X
        if g.args[0].grad != nil:
          g.args[0].grad = GraphValue[X](v:GraphValue[X](g.args[0].grad).v+t)
        else:
          g.args[0].grad = GraphValue[X](v:t)
      if gtConst notin g.args[1].tag:
        let t = GraphValue[R](g.grad).v.Y
        if g.args[1].grad != nil:
          g.args[1].grad = GraphValue[Y](v:GraphValue[Y](g.args[1].grad).v+t)
        else:
          g.args[1].grad = GraphValue[Y](v:t)))
proc `+`[X](x:GraphNode[X], y:SomeNumber):auto = x + newConst(y)
proc `+`[Y](x:SomeNumber, y:GraphNode[Y]):auto = newConst(x) + y

proc `*`[X,Y](x:GraphNode[X], y:GraphNode[Y]):auto =
  type R = type(GraphValue[X](x.val).v*GraphValue[Y](y.val).v)
  GraphNode[R](isop:true, str:"*", args: @[x.Graph,y],
    run: (proc(g:Graph) =
      echo "# Run: ",x.str," * ",y.str
      let v = GraphValue[R](v:GraphValue[X](g.args[0].val).v*GraphValue[Y](g.args[1].val).v)
      g.val = v
      g.str &= "(=" & $v.v & ")"),
    initGrad: (proc(g:Graph) = g.grad = GraphValue[R](v:1.R)),
    back: (proc(g:Graph) =
      echo "# Back: ",x.str," * ",y.str
      if gtConst notin g.args[0].tag:
        let t = GraphValue[R](g.grad).v*GraphValue[Y](g.args[1].val).v
        if g.args[0].grad != nil:
          g.args[0].grad = GraphValue[X](v:GraphValue[X](g.args[0].grad).v+t)
        else:
          g.args[0].grad = GraphValue[X](v:t)
      if gtConst notin g.args[1].tag:
        let t = GraphValue[X](g.args[0].val).v*GraphValue[R](g.grad).v
        if g.args[1].grad != nil:
          g.args[1].grad = GraphValue[X](v:GraphValue[X](g.args[1].grad).v+t)
        else:
          g.args[1].grad = GraphValue[Y](v:t)))
proc `*`[X](x:GraphNode[X], y:SomeNumber):auto = x * newConst(y)
proc `*`[Y](x:SomeNumber, y:GraphNode[Y]):auto = newConst(x) * y

when isMainModule:
  let
    x = newVar(2.0, "x")
    y = newVar(3.0, "y")
    z = x*(y+5.0)
    t = x*y*(z+1.0)*z
  echo "x: ",x
  echo "y: ",y
  echo "z: ",z
  echo "t: ",t
  let rt = t.eval
  echo "t: ",t
  echo "rt: ",rt
  echo "rz: ",z.eval
  echo "dtdz: ",t.evalGrad z
  echo "dtdx: ",t.evalGrad x
  echo "dtdy: ",t.evalGrad y
  echo "dzdx: ",z.evalGrad x
  echo "dzdy: ",z.evalGrad y
  let u = (x+t)*(t+z)
  echo "u: ",u
  echo "dudy: ",u.evalGrad y
