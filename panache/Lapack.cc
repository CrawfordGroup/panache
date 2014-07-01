#include <stdexcept> // for std::invalid_argument
#include "panache/Lapack.h" // for definition of int_t

#if FC_SYMBOL==1
#define F_DCOPY  dcopy
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
extern void F_DCOPY(int_t *length, double *x, int_t *inc_x,  double *y, int_t *inc_y);
extern double F_DDOT(int_t *n, double *x, int_t *inc_x, double *y, int_t *inc_y);
extern void F_DSCAL(int_t *n, double *alpha, double *vec, int_t *inc);
extern void F_DSWAP(int_t *length, double *x, int_t *inc_x, double *y, int_t *inc_y);

// Blas 2/3
extern void F_DGEMM(char*, char*, int_t*, int_t*, int_t*, double*, double*, int_t*, double*, int_t*, double*, double*, int_t*);
extern void F_DTRSM(char*, char*, char*, char*, int_t*, int_t*, double*, double*, int_t*, double*, int_t*);
extern void F_DSYMM(char*, char*, int_t*, int_t*, double*, double*, int_t*, double*, int_t*, double*, double*, int_t*);


// Lapack
extern int_t F_DGEQRF(int_t*, int_t*, double*, int_t*, double*, double*, int_t*,  int_t*);
extern int_t F_DORGQR(int_t*, int_t*, int_t*, double*, int_t*, double*, double*, int_t*,  int_t*);
extern int_t F_DPOTRF(char*, int_t*, double*, int_t*,  int_t*);
extern int_t F_DPOTRI(char*, int_t*, double*, int_t*,  int_t*);
extern int_t F_DSYEV(char*, char*, int_t*, double*, int_t*, double*, double*, int_t*,  int_t*);
} // end extern "C"



namespace panache {

// BLAS
void C_DCOPY(int_t length, double *x, int_t inc_x, double *y, int_t inc_y)
{
        ::F_DCOPY(&length, x, &inc_x, y, &inc_y);
}



double C_DDOT(int_t length, double *x, int_t inc_x, double *y, int_t inc_y)
{
    if(length == 0) return 0.0;

    return ::F_DDOT(&length, x, &inc_x, y, &inc_y);
}


void C_DSCAL(int_t length, double alpha, double *vec, int_t inc)
{
        ::F_DSCAL(&length, &alpha, vec, &inc);
}

void C_DSWAP(int_t length, double *x, int_t inc_x, double *y, int_t inc_y)
{
        ::F_DSWAP(&length, x, &inc_x, y, &inc_y);
}





// BLAS 2&3
void C_DGEMM(char transa, char transb, int_t m, int_t n, int_t k, 
             double alpha, double* a, int_t lda,
             double* b, int_t ldb, 
             double beta, double* c, int_t ldc)
{
    if(m == 0 || n == 0 || k == 0) return;
    ::F_DGEMM(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc);
}



void C_DTRSM(char side, char uplo, char transa, char diag,
             int_t m, int_t n, double alpha, double* a, int_t lda, double* b, int_t ldb)
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


void C_DSYMM(char side, char uplo, int_t m, int_t n, double alpha, double* a, int_t lda, double* b, int_t ldb, double beta, double* c, int_t ldc)
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
int_t C_DGEQRF(int_t m, int_t n, double* a, int_t lda, double* tau, double* work, int_t lwork)
{
    int_t info;
    ::F_DGEQRF(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}


int_t C_DORGQR(int_t m, int_t n, int_t k, double* a, int_t lda, double* tau, double* work, int_t lwork)
{
    int_t info;
    ::F_DORGQR(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

int_t C_DPOTRF(char uplo, int_t n, double* a, int_t lda)
{
    int_t info;
    ::F_DPOTRF(&uplo, &n, a, &lda, &info);
    return info;
}

int_t C_DPOTRI(char uplo, int_t n, double* a, int_t lda)
{
    int_t info;
    ::F_DPOTRI(&uplo, &n, a, &lda, &info);
    return info;
}

int_t C_DSYEV(char jobz, char uplo, int_t n, double* a, int_t lda, double* w, double* work, int_t lwork)
{
    int_t info;
    ::F_DSYEV(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
    return info;
}

} // close namespace panache
