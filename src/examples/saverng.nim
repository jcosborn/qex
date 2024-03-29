import qex

qexinit()

let
  lat = intSeqParam("lat", @[8,8,8,8])
  lo = lat.newLayout

var seed = 987654321u
var r = lo.newRNGField(RngMilc6, seed)
var r2 = lo.newRNGField(RngMilc6, seed+1)

proc check(x,y: Field) =
  for i in x:
    if x[i] != y[i]:
      echo "not equal ", i
      return
  echo "equal"

echo "should fail:"
check r, r2

const fileMd = "<?xml version=\"1.0\"?>\n<note>generated by QEX</note>\n"
const recordMd = "<?xml version=\"1.0\"?>\n<note>RNG field</note>\n"
var fn = "rng.lat"

var wr = r.l.newWriter(fn, fileMd)
wr.write(r, recordMd)
wr.close

var rd = r2.l.newReader(fn)
rd.read(r2)
rd.close

echo "should pass:"
check r, r2

qexfinalize()
