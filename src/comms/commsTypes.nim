type
  Comm* = ref object of RootObj

# globals

var defaultComm*: Comm
template getDefaultComm*(): Comm = defaultComm
template getComm*(): Comm = getDefaultComm()  # temporary alias
var myRank* = 0
var nRanks* = 1
