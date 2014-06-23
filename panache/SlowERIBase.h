#ifndef PANACHE_SLOWERIBASE_H
#define PANACHE_SLOWERIBASE_H

#define MAXFAC 100
#define EPS 1.0E-17


class SlowERIBase
{
private:
    double * df;
    double * fac;
    double ** bc;

    void calc_f(double *F, int n, double t);

    // Helper functions - mostly duplicates of stuff in libciomr
    double norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
                      double alpha1, const double* A);
    double* init_array(unsigned long int size);
    void free_array(double* array);


public:
SlowERIBase();
~SlowERIBase();

double eri(int l1, int m1, int n1, double alpha1,
           const double* A, int l2, int m2, int n2,
           double alpha2, const double* B, int l3, int m3,
           int n3, double alpha3, const double* C, int l4,
           int m4, int n4, double alpha4, const double* D,
           int norm_flag);

};

#endif

