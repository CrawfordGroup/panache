#include <stdexcept> // for std::invalid_argument
#include "Lapack.h" // for definition of int_t

#if FC_SYMBOL==1
#define F_DAXPY  daxpy
#define F_DCOPY  dcopy
#define F_DDOT   ddot
#define F_DROT   drot
#define F_DSCAL  dscal
#define F_DSWAP  dswap

#define F_DGEMM  dgemm
#define F_DGEMV  dgemv
#define F_DTRSM  dtrsm
#define F_DSYMM  dsymm
#define F_DTRMM  dtrmm

#define F_DGEQRF dgeqrf
#define F_DGESDD dgesdd
#define F_DGETRF dgetrf
#define F_DGETRI dgetri
#define F_DGETRS dgetrs
#define F_DORGQR dorgqr
#define F_DPOTRF dpotrf
#define F_DPOTRI dpotri
#define F_DSYEV  dsyev
#define F_DSYGV  dsygv


#elif FC_SYMBOL==2
#define F_DAXPY  daxpy_
#define F_DCOPY  dcopy_
#define F_DDOT   ddot_
#define F_DROT   drot_
#define F_DSCAL  dscal_
#define F_DSWAP  dswap_

#define F_DGEMM  dgemm_
#define F_DGEMV  dgemv_
#define F_DTRSM  dtrsm_
#define F_DSYMM  dsymm_
#define F_DTRMM  dtrmm_

#define F_DGEQRF dgeqrf_
#define F_DGESDD dgesdd_
#define F_DGETRF dgetrf_
#define F_DGETRI dgetri_
#define F_DGETRS dgetrs_
#define F_DORGQR dorgqr_
#define F_DPOTRF dpotrf_
#define F_DPOTRI dpotri_
#define F_DSYEV  dsyev_
#define F_DSYGV  dsygv_


#elif FC_SYMBOL==3
#define F_DAXPY  DAXPY
#define F_DCOPY  DCOPY
#define F_DDOT   DDOT
#define F_DROT   DROT
#define F_DSCAL  DSCAL
#define F_DSWAP  DSWAP

#define F_DGEMM  DGEMM
#define F_DGEMV  DGEMV
#define F_DTRSM  DTRSM
#define F_DSYMM  DSYMM
#define F_DTRMM  DTRMM

#define F_DGEQRF DGEQRF
#define F_DGESDD DGESDD
#define F_DGETRF DGETRF
#define F_DGETRI DGETRI
#define F_DGETRS DGETRS
#define F_DORGQR DORGQR
#define F_DPOTRF DPOTRF
#define F_DPOTRI DPOTRI
#define F_DSYEV  DSYEV
#define F_DSYGV  DSYGV


#elif FC_SYMBOL==4
#define F_DAXPY  DAXPY_
#define F_DCOPY  DCOPY_
#define F_DDOT   DDOT_
#define F_DROT   DROT_
#define F_DSCAL  DSCAL_
#define F_DSWAP  DSWAP_

#define F_DGEMM  DGEMM_
#define F_DGEMV  DGEMV_
#define F_DTRSM  DTRSM_
#define F_DSYMM  DSYMM_
#define F_DTRMM  DTRMM_

#define F_DGEQRF DGEQRF_
#define F_DGESDD DGESDD_
#define F_DGETRF DGETRF_
#define F_DGETRI DGETRI_
#define F_DGETRS DGETRS_
#define F_DORGQR DORGQR_
#define F_DPOTRF DPOTRF_
#define F_DPOTRI DPOTRI_
#define F_DSYEV  DSYEV_
#define F_DSYGV  DSYGV_


#endif



extern "C"
{

// Blas
extern void F_DAXPY(int_t *length, double *a, double *x, int_t *inc_x, double *y, int_t *inc_y);
extern void F_DCOPY(int_t *length, double *x, int_t *inc_x,  double *y, int_t *inc_y);
extern double F_DDOT(int_t *n, double *x, int_t *inc_x, double *y, int_t *inc_y);
extern void F_DROT(int_t *ntot, double *x, int_t *inc_x, double *y, int_t *inc_y, double *cotheta, double *sintheta);
extern void F_DSCAL(int_t *n, double *alpha, double *vec, int_t *inc);
extern void F_DSWAP(int_t *length, double *x, int_t *inc_x, double *y, int_t *inc_y);

// Blas 2/3
extern void F_DGEMM(char*, char*, int_t*, int_t*, int_t*, double*, double*, int_t*, double*, int_t*, double*, double*, int_t*);
extern void F_DGEMV(char*, int_t*, int_t*, double*, double*, int_t*, double*, int_t*, double*, double*, int_t*);
extern void F_DTRSM(char*, char*, char*, char*, int_t*, int_t*, double*, double*, int_t*, double*, int_t*);
extern void F_DSYMM(char*, char*, int_t*, int_t*, double*, double*, int_t*, double*, int_t*, double*, double*, int_t*);
extern void F_DTRMM(char*, char*, char*, char*, int_t*, int_t*, double*, double*, int_t*, double*, int_t*);


// Lapack
extern int_t F_DGEQRF(int_t*, int_t*, double*, int_t*, double*, double*, int_t*,  int_t*);
extern int_t F_DGESDD(char*, int_t*, int_t*, double*, int_t*, double*, double*, int_t*, double*, int_t*, double*, int_t*, int_t*,  int_t*);
extern int_t F_DGETRF(int_t*, int_t*, double*, int_t*, int_t*,  int_t*);
extern int_t F_DGETRI(int_t*, double*, int_t*, int_t*, double*, int_t*,  int_t*);
extern int_t F_DGETRS(char*, int_t*, int_t*, double*, int_t*, int_t*, double*, int_t*,  int_t*);
extern int_t F_DORGQR(int_t*, int_t*, int_t*, double*, int_t*, double*, double*, int_t*,  int_t*);
extern int_t F_DPOTRF(char*, int_t*, double*, int_t*,  int_t*);
extern int_t F_DPOTRI(char*, int_t*, double*, int_t*,  int_t*);
extern int_t F_DSYEV(char*, char*, int_t*, double*, int_t*, double*, double*, int_t*,  int_t*);
extern int_t F_DSYGV(int_t*, char*, char*, int_t*, double*, int_t*, double*, int_t*, double*, double*, int_t*,  int_t*);
} // end extern "C"



namespace panache {

// BLAS
void C_DAXPY(int_t length, double a, double *x, int_t inc_x, double *y, int_t inc_y)
{
    ::F_DAXPY(&length, &a, x, &inc_x, y, &inc_y);
}



void C_DCOPY(int_t length, double *x, int_t inc_x, double *y, int_t inc_y)
{
        ::F_DCOPY(&length, x, &inc_x, y, &inc_y);
}



double C_DDOT(int_t length, double *x, int_t inc_x, double *y, int_t inc_y)
{
    if(length == 0) return 0.0;

    return ::F_DDOT(&length, x, &inc_x, y, &inc_y);
}


void C_DROT(int_t length, double *x, int_t inc_x, 
            double *y, int_t inc_y, double costheta, double sintheta)
{
        ::F_DROT(&length, x, &inc_x, y, &inc_y, &costheta, &sintheta);
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



void C_DGEMV(char trans, int_t m, int_t n, 
             double alpha, double* a, int_t lda,
             double* x, int_t inc_x,
             double beta, double* y, int_t inc_y)
{
    if(m == 0 || n == 0) return;
    if (trans == 'N' || trans == 'n') trans = 'T';
    else if (trans == 'T' || trans == 't') trans = 'N';
    else throw std::invalid_argument("C_DGEMV trans argument is invalid.");
    ::F_DGEMV(&trans, &n, &m, &alpha, a, &lda, x, &inc_x, &beta, y, &inc_y);
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


void C_DTRMM(char side, char uplo, char transa, char diag, int_t m, int_t n, double alpha, double* a, int_t lda, double* b, int_t ldb)
{
    if(m == 0 || n == 0) return;
    if (uplo == 'U' || uplo == 'u') uplo = 'L';
    else if (uplo == 'L' || uplo == 'l') uplo = 'U';
    else throw std::invalid_argument("C_DTRMM uplo argument is invalid.");
    if (side == 'L' || side == 'L') side = 'R';
    else if (side == 'R' || side == 'r') side = 'L';
    else throw std::invalid_argument("C_DTRMM side argument is invalid.");
    ::F_DTRMM(&side, &uplo, &transa, &diag, &n, &m, &alpha, a, &lda, b, &ldb);
}



// Lapack
int_t C_DGEQRF(int_t m, int_t n, double* a, int_t lda, double* tau, double* work, int_t lwork)
{
    int_t info;
    ::F_DGEQRF(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}


int_t C_DGESDD(char jobz, int_t m, int_t n, double* a, int_t lda, 
             double* s, double* u, int_t ldu,
             double* vt, int_t ldvt, double* work, int_t lwork, int_t* iwork)
{
    int_t info;
    ::F_DGESDD(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);
    return info;
}


int_t C_DGETRF(int_t m, int_t n, double* a, int_t lda, int_t* ipiv)
{
    int_t info;
    ::F_DGETRF(&m, &n, a, &lda, ipiv, &info);
    return info;
}

int_t C_DGETRI(int_t n, double* a, int_t lda, int_t* ipiv, double* work, int_t lwork)
{
    int_t info;
    ::F_DGETRI(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

int_t C_DGETRS(char trans, int_t n, int_t nrhs, double* a, int_t lda, int_t* ipiv, double* b, int_t ldb)
{
    int_t info;
    ::F_DGETRS(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
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

int_t C_DSYGV(int_t itype, char jobz, char uplo, int_t n,
            double* a, int_t lda,
            double* b, int_t ldb,
            double* w, double* work, int_t lwork)
{
    int_t info;
    ::F_DSYGV(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info);
    return info;
}

} // close namespace panache
