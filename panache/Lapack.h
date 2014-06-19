#ifndef PANACHE_LAPACK_H
#define PANACHE_LAPACK_H

#include "int_t.h"


namespace panache {


// BLAS
void C_DAXPY(lapack_int_t length, double a, double *x, lapack_int_t inc_x, double *y, lapack_int_t inc_y);

void C_DCOPY(lapack_int_t length, double *x, lapack_int_t inc_x, double *y, lapack_int_t inc_y);

double C_DDOT(lapack_int_t length, double *X, lapack_int_t inc_x, double *Y, lapack_int_t inc_y);

void C_DROT(lapack_int_t length, double *x, lapack_int_t inc_x, 
            double *y, lapack_int_t inc_y, double costheta, double sintheta);

void C_DSCAL(lapack_int_t length, double alpha, double *vec, lapack_int_t inc);

void C_DSWAP(lapack_int_t length, double *x, lapack_int_t inc_x, double *y, lapack_int_t inc_y);


// BLAS 2/3
void C_DGEMM(char transa, char transb, lapack_int_t m, lapack_int_t n, lapack_int_t k, 
             double alpha, double* a, lapack_int_t lda,
             double* b, lapack_int_t ldb, 
             double beta, double* c, lapack_int_t ldc);



void C_DGEMV(char trans, lapack_int_t m, lapack_int_t n, 
             double alpha, double* a, lapack_int_t lda,
             double* x, lapack_int_t inc_x,
             double beta, double* y, lapack_int_t inc_y);

void C_DTRSM(char side, char uplo, char transa, char diag,
             lapack_int_t m, lapack_int_t n, double alpha, double* a, lapack_int_t lda, double* b, lapack_int_t ldb);


void C_DSYMM(char side, char uplo, lapack_int_t m, lapack_int_t n, double alpha,
             double* a, lapack_int_t lda, double* b, lapack_int_t ldb,
             double beta, double* c, lapack_int_t ldc);

void C_DTRMM(char side, char uplo, char transa, char diag,
             lapack_int_t m, lapack_int_t n, double alpha, double* a, lapack_int_t lda, double * b, lapack_int_t ldb);



// Lapack
lapack_int_t C_DGEQRF(lapack_int_t m, lapack_int_t n, double* a, lapack_int_t lda, double* tau, double* work, lapack_int_t lwork);

lapack_int_t C_DGESDD(char jobz, lapack_int_t m, lapack_int_t n, double* a, lapack_int_t lda, 
             double* s, double* u, lapack_int_t ldu,
             double* vt, lapack_int_t ldvt, double* work, lapack_int_t lwork, lapack_int_t * iwork);


lapack_int_t C_DGETRF(lapack_int_t m, lapack_int_t n, double* a, lapack_int_t lda, lapack_int_t * ipiv);

lapack_int_t C_DGETRI(lapack_int_t n, double* a, lapack_int_t lda, lapack_int_t * ipiv, double* work, lapack_int_t lwork);

lapack_int_t C_DGETRS(char trans, lapack_int_t n, lapack_int_t nrhs, double* a, lapack_int_t lda, lapack_int_t * ipiv, double* b, lapack_int_t ldb);

lapack_int_t C_DORGQR(lapack_int_t m, lapack_int_t n, lapack_int_t k, double* a, lapack_int_t lda, double* tau, double* work, lapack_int_t lwork);

lapack_int_t C_DPOTRF(char uplo, lapack_int_t n, double* a, lapack_int_t lda);

lapack_int_t C_DPOTRI(char uplo, lapack_int_t n, double* a, lapack_int_t lda);




lapack_int_t C_DSYEV(char jobz, char uplo, lapack_int_t n, double* a, lapack_int_t lda, double* w, double* work, lapack_int_t lwork);

lapack_int_t C_DSYGV(lapack_int_t itype, char jobz, char uplo, lapack_int_t n,
            double* a, lapack_int_t lda,
            double* b, lapack_int_t ldb,
            double* w, double* work, lapack_int_t lwork);
}

#endif
