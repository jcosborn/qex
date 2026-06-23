import ../[core, scalar, functional]

let grt = initGraphRuntime()

block:
  let x = grt.toGvalue(3.0)
  let y = grt.toGvalue(2.0)
  let v = grt.localScalar()
  let z = apply(lambda(v, v + v), x * y)
  echo "## z (before eval)"
  echo z.treeRepr
  z.eval
  echo "## z (after eval)"
  echo z.treeRepr
  echo "z = ", z
  let dzdx = z.grad x
  echo "## dzdx (before eval)"
  echo dzdx.treeRepr
  dzdx.eval
  echo "## dzdx (after eval)"
  echo dzdx.treeRepr
  echo "dzdx = ", dzdx

block:
  let f = lambdaParam(grt.localScalar(), grt.localScalar())
  let u = grt.localScalar()
  let hof = lambda(f, lambda(u, Gscalar(apply(f, u)) + 1.0))
  let v = grt.localScalar()
  let a = grt.toGvalue(2.0)
  let g = lambda(v, a * v)
  let z = apply(apply(hof, g), 3.0)
  echo "## z (before eval)"
  echo z.treeRepr
  z.eval
  echo "## z (after eval)"
  echo z.treeRepr
  echo "z = ", z
  let dzda = z.grad a
  echo "## dzda (before eval)"
  echo dzda.treeRepr
  dzda.eval
  echo "## dzda (after eval)"
  echo dzda.treeRepr
  echo "dzda = ", dzda

block:
  # Short recursion via the Y combinator.
  let protoArg = grt.localScalar()
  let protoRet = grt.localScalar()
  let fnProto = lambda(protoArg, protoRet)
  let x = lambdaParam(fnProto, fnProto)
  let f = lambdaParam(fnProto, fnProto)
  let Y = lambda(f, apply(lambda(x, apply(f, apply(x, x))), lambda(x, apply(f, apply(x, x)))))

  let rf = lambdaParam(grt.localScalar(), grt.localScalar())
  let u = grt.localScalar()
  let v = grt.localScalar()
  let base = grt.toGvalue(1.0)
  let step = grt.toGvalue(1.0)
  let F = lambda(rf, lambda(u,
    cond(equal(u, 0.0), base,
      Gscalar(apply(lambda(v, Gscalar(apply(rf, v)) + step), u - 1.0)))))

  let z = apply(apply(Y, F), 3.0)
  echo "## Y recursion z (before eval)"
  echo z.treeRepr
  z.eval
  echo "## Y recursion z (after eval)"
  echo z.treeRepr
  echo "Y recursion z = ", z
  let dzdstep = z.grad step
  echo "## Y recursion dzdstep (before eval)"
  echo dzdstep.treeRepr
  dzdstep.eval
  echo "## Y recursion dzdstep (after eval)"
  echo dzdstep.treeRepr
  echo "Y recursion dzdstep = ", dzdstep
  let dzdbase = z.grad base
  echo "## Y recursion dzdbase (before eval)"
  echo dzdbase.treeRepr
  dzdbase.eval
  echo "## Y recursion dzdbase (after eval)"
  echo dzdbase.treeRepr
  echo "Y recursion dzdbase = ", dzdbase
