#ifndef PANACHE_LAPACK_H
#define PANACHE_LAPACK_H

#include <cstdint>
#include <cstdio> // will not be needed once schmidt and invert_matrix prototypes are moved
#ifdef USE_64LAPACK
  #define int_t int64_t
#else
  #define int_t int32_t
#endif

namespace panache {

//! \todo MOVE THESE
void schmidt(double **A, int rows, int cols, FILE *outfile);
double invert_matrix(double **a, double **y, int N, FILE *outfile);



// BLAS
void C_DAXPY(int_t length, double a, double *x, int_t inc_x, double *y, int_t inc_y);

void C_DCOPY(int_t length, double *x, int_t inc_x, double *y, int_t inc_y);

double C_DDOT(int_t length, double *X, int_t inc_x, double *Y, int_t inc_y);

void C_DROT(int_t length, double *x, int_t inc_x, 
            double *y, int_t inc_y, double costheta, double sintheta);

void C_DSCAL(int_t length, double alpha, double *vec, int_t inc);

void C_DSWAP(int_t length, double *x, int_t inc_x, double *y, int_t inc_y);


// BLAS 2/3
void C_DGEMM(char transa, char transb, int_t m, int_t n, int_t k, 
             double alpha, double* a, int_t lda,
             double* b, int_t ldb, 
             double beta, double* c, int_t ldc);



void C_DGEMV(char trans, int_t m, int_t n, 
             double alpha, double* a, int_t lda,
             double* x, int_t inc_x,
             double beta, double* y, int_t inc_y);

void C_DTRSM(char side, char uplo, char transa, char diag,
             int_t m, int_t n, double alpha, double* a, int_t lda, double* b, int_t ldb);




// Lapack
int_t C_DGEQRF(int_t m, int_t n, double* a, int_t lda, double* tau, double* work, int_t lwork);

int_t C_DGESDD(char jobz, int_t m, int_t n, double* a, int_t lda, 
             double* s, double* u, int_t ldu,
             double* vt, int_t ldvt, double* work, int_t lwork, int_t * iwork);


int_t C_DGETRF(int_t m, int_t n, double* a, int_t lda, int_t * ipiv);

int_t C_DGETRI(int_t n, double* a, int_t lda, int_t * ipiv, double* work, int_t lwork);

int_t C_DGETRS(char trans, int_t n, int_t nrhs, double* a, int_t lda, int_t * ipiv, double* b, int_t ldb);

int_t C_DORGQR(int_t m, int_t n, int_t k, double* a, int_t lda, double* tau, double* work, int_t lwork);

int_t C_DPOTRF(char uplo, int_t n, double* a, int_t lda);

int_t C_DPOTRI(char uplo, int_t n, double* a, int_t lda);




int_t C_DSYEV(char jobz, char uplo, int_t n, double* a, int_t lda, double* w, double* work, int_t lwork);

int_t C_DSYGV(int_t itype, char jobz, char uplo, int_t n,
            double* a, int_t lda,
            double* b, int_t ldb,
            double* w, double* work, int_t lwork);
}

#endif
