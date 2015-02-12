#include <stdexcept> // for std::invalid_argument
#include "panache/Lapack.h" // for definition of lapack_int_t

#if FC_SYMBOL==1
#define F_DCOPY  dcopy
#define F_DAXPY  daxpy
#define F_DDOT   ddot
#define F_DSCAL  dscal
#define F_DSWAP  dswap

#define F_DGEMM  dgemm
#define F_DTRSM  dtrsm
#define F_DSYMM  dsymm

#define F_DGEQRF dgeqrf
#define F_DORGQR dorgqr
#define F_DPOTRF dpotrf
#define F_DPOTRI dpotri
#define F_DSYEV  dsyev


#elif FC_SYMBOL==2
#define F_DCOPY  dcopy_
#define F_DAXPY  daxpy_
#define F_DDOT   ddot_
#define F_DSCAL  dscal_
#define F_DSWAP  dswap_

#define F_DGEMM  dgemm_
#define F_DTRSM  dtrsm_
#define F_DSYMM  dsymm_

#define F_DGEQRF dgeqrf_
#define F_DORGQR dorgqr_
#define F_DPOTRF dpotrf_
#define F_DPOTRI dpotri_
#define F_DSYEV  dsyev_


#elif FC_SYMBOL==3
#define F_DCOPY  DCOPY
#define F_DAXPY  DAXPY
#define F_DDOT   DDOT
#define F_DSCAL  DSCAL
#define F_DSWAP  DSWAP

#define F_DGEMM  DGEMM
#define F_DTRSM  DTRSM
#define F_DSYMM  DSYMM

#define F_DGEQRF DGEQRF
#define F_DORGQR DORGQR
#define F_DPOTRF DPOTRF
#define F_DPOTRI DPOTRI
#define F_DSYEV  DSYEV


#elif FC_SYMBOL==4
#define F_DCOPY  DCOPY_
#define F_DAXPY  DAXPY_
#define F_DDOT   DDOT_
#define F_DSCAL  DSCAL_
#define F_DSWAP  DSWAP_

#define F_DGEMM  DGEMM_
#define F_DTRSM  DTRSM_
#define F_DSYMM  DSYMM_

#define F_DGEQRF DGEQRF_
#define F_DORGQR DORGQR_
#define F_DPOTRF DPOTRF_
#define F_DPOTRI DPOTRI_
#define F_DSYEV  DSYEV_


#endif



extern "C"
{

// Blas
extern void F_DCOPY(lapack_int_t *length, double *x, lapack_int_t *inc_x,  double *y, lapack_int_t *inc_y);
extern void F_DAXPY(lapack_int_t *length, double *a, double *x, lapack_int_t *inc_x,  double *y, lapack_int_t *inc_y);
extern double F_DDOT(lapack_int_t *n, double *x, lapack_int_t *inc_x, double *y, lapack_int_t *inc_y);
extern void F_DSCAL(lapack_int_t *n, double *alpha, double *vec, lapack_int_t *inc);
extern void F_DSWAP(lapack_int_t *length, double *x, lapack_int_t *inc_x, double *y, lapack_int_t *inc_y);

// Blas 2/3
extern void F_DGEMM(char*, char*, lapack_int_t*, lapack_int_t*, lapack_int_t*, double*, double*, lapack_int_t*, double*, lapack_int_t*, double*, double*, lapack_int_t*);
extern void F_DTRSM(char*, char*, char*, char*, lapack_int_t*, lapack_int_t*, double*, double*, lapack_int_t*, double*, lapack_int_t*);
extern void F_DSYMM(char*, char*, lapack_int_t*, lapack_int_t*, double*, double*, lapack_int_t*, double*, lapack_int_t*, double*, double*, lapack_int_t*);


// Lapack
extern lapack_int_t F_DGEQRF(lapack_int_t*, lapack_int_t*, double*, lapack_int_t*, double*, double*, lapack_int_t*,  lapack_int_t*);
extern lapack_int_t F_DORGQR(lapack_int_t*, lapack_int_t*, lapack_int_t*, double*, lapack_int_t*, double*, double*, lapack_int_t*,  lapack_int_t*);
extern lapack_int_t F_DPOTRF(char*, lapack_int_t*, double*, lapack_int_t*,  lapack_int_t*);
extern lapack_int_t F_DPOTRI(char*, lapack_int_t*, double*, lapack_int_t*,  lapack_int_t*);
extern lapack_int_t F_DSYEV(char*, char*, lapack_int_t*, double*, lapack_int_t*, double*, double*, lapack_int_t*,  lapack_int_t*);
} // end extern "C"



namespace panache {

// BLAS
void C_DCOPY(lapack_int_t length, double *x, lapack_int_t inc_x, double *y, lapack_int_t inc_y)
{
        ::F_DCOPY(&length, x, &inc_x, y, &inc_y);
}


void C_DAXPY(lapack_int_t length, double a, double *x, lapack_int_t inc_x, double *y, lapack_int_t inc_y)
{
        ::F_DAXPY(&length, &a, x, &inc_x, y, &inc_y);
}


double C_DDOT(lapack_int_t length, double *x, lapack_int_t inc_x, double *y, lapack_int_t inc_y)
{
    if(length == 0) return 0.0;

    return ::F_DDOT(&length, x, &inc_x, y, &inc_y);
}


void C_DSCAL(lapack_int_t length, double alpha, double *vec, lapack_int_t inc)
{
        ::F_DSCAL(&length, &alpha, vec, &inc);
}

void C_DSWAP(lapack_int_t length, double *x, lapack_int_t inc_x, double *y, lapack_int_t inc_y)
{
        ::F_DSWAP(&length, x, &inc_x, y, &inc_y);
}





// BLAS 2&3
void C_DGEMM(char transa, char transb, lapack_int_t m, lapack_int_t n, lapack_int_t k, 
             double alpha, double* a, lapack_int_t lda,
             double* b, lapack_int_t ldb, 
             double beta, double* c, lapack_int_t ldc)
{
    if(m == 0 || n == 0 || k == 0) return;
    ::F_DGEMM(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc);
}



void C_DTRSM(char side, char uplo, char transa, char diag,
             lapack_int_t m, lapack_int_t n, double alpha, double* a, lapack_int_t lda, double* b, lapack_int_t ldb)
{
    if(m == 0 || n == 0) return;
    if (uplo == 'U' || uplo == 'u') uplo = 'L';
    else if (uplo == 'L' || uplo == 'l') uplo = 'U';
    else throw std::invalid_argument("C_DTRSM uplo argument is invalid.");
    if (side == 'L' || side == 'L') side = 'R';
    else if (side == 'R' || side == 'r') side = 'L';
    else throw std::invalid_argument("C_DTRSM side argument is invalid.");
    ::F_DTRSM(&side, &uplo, &transa, &diag, &n, &m, &alpha, a, &lda, b, &ldb);
}


void C_DSYMM(char side, char uplo, lapack_int_t m, lapack_int_t n, double alpha, double* a, lapack_int_t lda, double* b, lapack_int_t ldb, double beta, double* c, lapack_int_t ldc)
{
    if(m == 0 || n == 0) return;
    if (uplo == 'U' || uplo == 'u') uplo = 'L';
    else if (uplo == 'L' || uplo == 'l') uplo = 'U';
    else throw std::invalid_argument("C_DSYMM uplo argument is invalid.");
    if (side == 'L' || side == 'L') side = 'R';
    else if (side == 'R' || side == 'r') side = 'L';
    else throw std::invalid_argument("C_DSYMM side argument is invalid.");
    ::F_DSYMM(&side, &uplo, &n, &m, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}


// Lapack
lapack_int_t C_DGEQRF(lapack_int_t m, lapack_int_t n, double* a, lapack_int_t lda, double* tau, double* work, lapack_int_t lwork)
{
    lapack_int_t info;
    ::F_DGEQRF(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}


lapack_int_t C_DORGQR(lapack_int_t m, lapack_int_t n, lapack_int_t k, double* a, lapack_int_t lda, double* tau, double* work, lapack_int_t lwork)
{
    lapack_int_t info;
    ::F_DORGQR(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

lapack_int_t C_DPOTRF(char uplo, lapack_int_t n, double* a, lapack_int_t lda)
{
    lapack_int_t info;
    ::F_DPOTRF(&uplo, &n, a, &lda, &info);
    return info;
}

lapack_int_t C_DPOTRI(char uplo, lapack_int_t n, double* a, lapack_int_t lda)
{
    lapack_int_t info;
    ::F_DPOTRI(&uplo, &n, a, &lda, &info);
    return info;
}

lapack_int_t C_DSYEV(char jobz, char uplo, lapack_int_t n, double* a, lapack_int_t lda, double* w, double* work, lapack_int_t lwork)
{
    lapack_int_t info;
    ::F_DSYEV(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
    return info;
}

} // close namespace panache
