const EPS = 1e-15

# I0(x) = sum_{k=0 to infinity} ( (x/2)^(2k) ) / (k!)^2
proc besselI0*(x: float): float =
  let halfX = x / 2.0
  var
    term = 1.0   # term for k = 0
    sumVal = 1.0 # sum starts with k = 0
    k = 1
  while abs(term) > EPS * abs(sumVal):
    # term *= ( (x/2)^2 ) / (k^2 )
    term *= (halfX * halfX) / (float(k) * float(k))
    sumVal += term
    inc k
  sumVal

# I1(x) = sum_{k=0 to infinity} ( (x/2)^(2k+1) ) / ( k! * (k+1)! )
proc besselI1*(x: float): float =
  let halfX = x / 2.0
  var
    term = halfX # term for k = 0
    sumVal = halfX
    k = 1
  while abs(term) > EPS * abs(sumVal):
    # term *= ( (x/2)^2 ) / [ k * (k+1) ]
    term *= (halfX * halfX) / (float(k) * float(k + 1))
    sumVal += term
    inc k
  sumVal

# In(x) = sum_{k=0 to infinity} ( (x/2)^(2k+n) ) / (k! * (n+k)!)
proc besselIn*(n_in: int, x: float): float =
  let n = abs(n_in)
  let halfX = x / 2.0
  if x == 0.0:
    if n == 0: return 1.0
    else: return 0.0
  var term = 1.0
  for i in 1..n:
    term *= halfX / float(i)
  var sumVal = term
  var k = 1
  while abs(term) > EPS * abs(sumVal):
    # term *= (x/2)^2 / (k * (n+k))
    term *= (halfX * halfX) / (float(k) * float(n + k))
    sumVal += term
    inc k
  sumVal
