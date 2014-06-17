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

double eri(unsigned int l1, unsigned int m1, unsigned int n1, double alpha1,
           const double* A, unsigned int l2, unsigned int m2, unsigned int n2,
           double alpha2, const double* B, unsigned int l3, unsigned int m3,
           unsigned int n3, double alpha3, const double* C, unsigned int l4,
           unsigned int m4, unsigned int n4, double alpha4, const double* D,
           int norm_flag);

};

#endif

