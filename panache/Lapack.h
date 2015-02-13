#ifndef PANACHE_LAPACK_H
#define PANACHE_LAPACK_H

#include <cstdint>

/*! \def lapack_int_t
 *
 *  \brief Is set to either int64_t or int32_t, depending on the PANACHE_USE_64LAPACK option.
 *         See \ref compiling_lapack64
 */
#ifdef PANACHE_LAPACK64
  #define lapack_int_t int64_t
#else
  #define lapack_int_t int32_t
#endif


namespace panache {


// BLAS

void C_DCOPY(lapack_int_t length, double *x, lapack_int_t inc_x, double *y, lapack_int_t inc_y);

double C_DDOT(lapack_int_t length, double *X, lapack_int_t inc_x, double *Y, lapack_int_t inc_y);

void C_DSCAL(lapack_int_t length, double alpha, double *vec, lapack_int_t inc);

void C_DSWAP(lapack_int_t length, double *x, lapack_int_t inc_x, double *y, lapack_int_t inc_y);

void C_DAXPY(lapack_int_t length, double a, double *x, lapack_int_t inc_x, double *y, lapack_int_t inc_y);


// BLAS 2/3
void C_DGEMM(char transa, char transb, lapack_int_t m, lapack_int_t n, lapack_int_t k, 
             double alpha, double* a, lapack_int_t lda,
             double* b, lapack_int_t ldb, 
             double beta, double* c, lapack_int_t ldc);



void C_DTRSM(char side, char uplo, char transa, char diag,
             lapack_int_t m, lapack_int_t n, double alpha, double* a, lapack_int_t lda, double* b, lapack_int_t ldb);


void C_DSYMM(char side, char uplo, lapack_int_t m, lapack_int_t n, double alpha,
             double* a, lapack_int_t lda, double* b, lapack_int_t ldb,
             double beta, double* c, lapack_int_t ldc);


// Lapack
lapack_int_t C_DGEQRF(lapack_int_t m, lapack_int_t n, double* a, lapack_int_t lda, double* tau, double* work, lapack_int_t lwork);

lapack_int_t C_DORGQR(lapack_int_t m, lapack_int_t n, lapack_int_t k, double* a, lapack_int_t lda, double* tau, double* work, lapack_int_t lwork);

lapack_int_t C_DPOTRF(char uplo, lapack_int_t n, double* a, lapack_int_t lda);

lapack_int_t C_DPOTRI(char uplo, lapack_int_t n, double* a, lapack_int_t lda);




lapack_int_t C_DSYEV(char jobz, char uplo, lapack_int_t n, double* a, lapack_int_t lda, double* w, double* work, lapack_int_t lwork);

} // close namespace panache

#endif
