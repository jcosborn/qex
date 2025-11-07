typedef int integer;
//typedef long long integer;
typedef double doublereal;

#define U(x) x ## _
#define L(r,f,a) \
  r f a; \
  r U(f) a

L(void, cgemm, (const char *transa, const char *transb,
		const int *m, const int *n, const int *k,
		const float *alpha, const float *a,
		const int *lda, const float *b, const int *ldb,
		const float *beta, float *c, const int *ldc));

L(void, dgemm, (const char *transa, const char *transb,
		const int *m, const int *n, const int *k,
		const double *alpha, const double *a,
		const int *lda, const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc));

L(void, dsterf, (integer *n, doublereal *d, doublereal *e, integer *info));
L(void, dgetrf, (integer *m, integer *n, doublereal *a, integer * lda,
		 integer *ipiv, integer *info));
L(void, dbdsqr, (char *uplo, integer *n, integer *ncvt, integer *nru,
		 integer *ncc, doublereal *d, doublereal *e, doublereal *vt,
		 integer *ldvt, doublereal *u, integer *ldu, doublereal *c,
		 integer *ldc, doublereal *work, integer *info));
