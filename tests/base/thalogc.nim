# thalo under the GC stress and borrow check of thalogc.nims.  makeGatherMap
# lends its send lists to scatter only as raw pointers.  If they are not kept
# referenced, the check aborts at scatter; without it, a build whose compiler
# keeps no stack copy sends freed memory, and "map indices" finds indices out
# of range.
include thalo
