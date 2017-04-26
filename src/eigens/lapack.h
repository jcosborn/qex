typedef int integer;
//typedef long long integer;
typedef double doublereal;

#define U(x) x ## _
#define L(r,f,a) \
  r f a; \
  r U(f) a

L(void, dsterf, (integer *n, doublereal *d, doublereal *e, integer *info));
L(void, dgetrf, (integer *m, integer *n, doublereal *a, integer * lda,
		 integer *ipiv, integer *info));
L(void, dbdsqr, (char *uplo, integer *n, integer *ncvt, integer *nru,
		 integer *ncc, doublereal *d, doublereal *e, doublereal *vt,
		 integer *ldvt, doublereal *u, integer *ldu, doublereal *c,
		 integer *ldc, doublereal *work, integer *info));
