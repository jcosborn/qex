# flops counts
# all names start with 'f' for flops

# Complex number (fc...)
template fcadd*:int = 2
template fcmul*:int = 6
template fcredot*:int = 3

# Complex matrix (fcm...)
template fcmadd*(n:int):int = n*n*fcadd
template fcmmul*(n:int):int = n*n*(n*fcmul + (n-1)*fcadd)
template fcmredot*(n:int):int = n*n*fcredot + (n*n-1)

# single plaquette: redot(A*B,C*D)
template fplaq*(n:int):int = 2*fcmmul(n) + fcmredot(n)
