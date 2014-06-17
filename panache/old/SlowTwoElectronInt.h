#ifndef PANACHE_SLOWTWOELECTRONINT_H
#define PANACHE_SLOWTWOELECTRONINT_H

#include "TwoBodyAOInt.h"

namespace panache
{


/*! \ingroup MINTS
 *  \class ERI
 *  \brief Capable of computing two-electron repulsion integrals.
 */
class SlowTwoElectronInt : public TwoBodyAOInt
{
private:
    double eri(unsigned int l1, unsigned int m1, unsigned int n1, double alpha1,
               const double* A, unsigned int l2, unsigned int m2, unsigned int n2,
               double alpha2, const double* B, unsigned int l3, unsigned int m3,
               unsigned int n3, double alpha3, const double* C, unsigned int l4,
               unsigned int m4, unsigned int n4, double alpha4, const double* D,
               int norm_flag);

    void calc_f(double *, int, double);

    double norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
                      double alpha1, const double* A);

    void free_array(double * arr)
    {
        delete [] arr;
    }

    double *df;
    double *fac;
    double **bc;



protected:
    //! Maximum cartesian class size.
    int max_cart_;

    //! Computes the ERIs between four shells.
    size_t compute_quartet(int, int, int, int);

    //! Original shell index requested
    int osh1_, osh2_, osh3_, osh4_;


public:
    //! Constructor
    SlowTwoElectronInt(const std::shared_ptr<BasisSet> bs1,
                       const std::shared_ptr<BasisSet> bs2,
                       const std::shared_ptr<BasisSet> bs3,
                       const std::shared_ptr<BasisSet> bs4);

    virtual ~SlowTwoElectronInt();

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(const AOShellCombinationsIterator&);

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    virtual size_t compute_shell(int, int, int, int);

};

}

#endif //PANACHE_TWOELECTRONINT_H


