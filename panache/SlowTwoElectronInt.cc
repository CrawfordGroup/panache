#include <cstring> // memset, etc
#include "SlowERIBase.h"
#include "SlowTwoElectronInt.h"
#include "BasisSet.h"
#include "BasisFunctionMacros.h"
#include "AOShellCombinationsIterator.h"
#include "CartesianIter.h"

namespace panache
{


SlowTwoElectronInt::SlowTwoElectronInt(const SharedBasisSet bs1,
                                       const SharedBasisSet bs2,
                                       const SharedBasisSet bs3,
                                       const SharedBasisSet bs4)
    : TwoBodyAOInt(bs1,bs2,bs3,bs4)
{
    // Note - there is no permutation, etc, so the bs#_ is the same as original_bs#_
    bs1_ = bs1;
    bs2_ = bs2;
    bs3_ = bs3;
    bs4_ = bs4;

    size_t size = INT_NCART(basis1()->max_am()) * INT_NCART(basis2()->max_am()) *
                  INT_NCART(basis3()->max_am()) * INT_NCART(basis4()->max_am());

    // Used in pure_transform
    tformbuf_ = new double[size];
    memset(tformbuf_, 0, sizeof(double)*size);

    target_ = new double[size];
    memset(target_, 0, sizeof(double)*size);

    source_ = new double[size];
    memset(source_, 0, sizeof(double)*size);
}

SlowTwoElectronInt::~SlowTwoElectronInt()
{
    delete[] tformbuf_;
    delete[] target_;
    delete[] source_;
}

size_t SlowTwoElectronInt::compute_shell(const AOShellCombinationsIterator& shellIter)
{
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}

size_t SlowTwoElectronInt::compute_shell(int sh1, int sh2, int sh3, int sh4)
{
    int n1, n2, n3, n4;

    if (force_cartesian_)
    {
        n1 = original_bs1_->shell(sh1).ncartesian();
        n2 = original_bs2_->shell(sh2).ncartesian();
        n3 = original_bs3_->shell(sh3).ncartesian();
        n4 = original_bs4_->shell(sh4).ncartesian();
    }
    else
    {
        n1 = original_bs1_->shell(sh1).nfunction();
        n2 = original_bs2_->shell(sh2).nfunction();
        n3 = original_bs3_->shell(sh3).nfunction();
        n4 = original_bs4_->shell(sh4).nfunction();
    }
    curr_buff_size_ = n1 * n2 * n3 * n4;

    size_t ncomputed = compute_quartet(sh1, sh2, sh3, sh4);

    if (ncomputed)
    {
        // Only do the following if we did any work.
        // copy the integrals to the target_
        std::copy(source_, source_ + curr_buff_size_, target_);
    }

    return ncomputed;
}

size_t SlowTwoElectronInt::compute_quartet(int sh1, int sh2, int sh3, int sh4)
{
    const GaussianShell& s1 = original_bs1_->shell(sh1);
    const GaussianShell& s2 = original_bs2_->shell(sh2);
    const GaussianShell& s3 = original_bs3_->shell(sh3);
    const GaussianShell& s4 = original_bs4_->shell(sh4);

    int am1 = s1.am();
    int am2 = s2.am();
    int am3 = s3.am();
    int am4 = s4.am();

    double A[3], B[3], C[3], D[3];

    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];
    C[0] = s3.center()[0];
    C[1] = s3.center()[1];
    C[2] = s3.center()[2];
    D[0] = s4.center()[0];
    D[1] = s4.center()[1];
    D[2] = s4.center()[2];


    CartesianIter cit1(am1);
    CartesianIter cit2(am2);
    CartesianIter cit3(am3);
    CartesianIter cit4(am4);

    int ncart1 = INT_NCART(am1);
    int ncart2 = INT_NCART(am2);
    int ncart3 = INT_NCART(am3);
    int ncart4 = INT_NCART(am4);
    size_t size = ncart1 * ncart2 * ncart3 * ncart4;


    std::array<int, 3> exp1, exp2, exp3, exp4;

    int ijkl = 0;

    cit1.start();
    for(int i = 0; i < ncart1; i++)
    {
        exp1 = cit1.exponents();
        cit2.start();

        for(int j = 0; j < ncart2; j++)
        {
            exp2 = cit2.exponents();
            cit3.start();

            for(int k = 0; k < ncart3; k++)
            {
                exp3 = cit3.exponents();
                cit4.start();

                for(int l = 0; l < ncart4; l++)
                {
                    exp4 = cit4.exponents();

                    source_[ijkl] = 0;

                    // Loop over primitives
                    for(int a = 0; a < s1.nprimitive(); a++)
                        for(int b = 0; b < s2.nprimitive(); b++)
                            for(int c = 0; c < s3.nprimitive(); c++)
                                for(int d = 0; d < s4.nprimitive(); d++)
                                    source_[ijkl] += sloweri_.eri(
                                                         exp1[0], exp1[1], exp1[2], s1.exp(a), A,
                                                         exp2[0], exp2[1], exp2[2], s2.exp(b), B,
                                                         exp3[0], exp3[1], exp3[2], s3.exp(c), C,
                                                         exp4[0], exp4[1], exp4[2], s4.exp(d), D,
                                                         0)*s1.coef(a)*s2.coef(b)*s3.coef(c)*s4.coef(d);
                    ijkl++;
                    cit4.next();
                }

                cit3.next();

            }

            cit2.next();

        }

        cit1.next();
    }

    // Transform the integrals into pure angular momentum
    if (!force_cartesian_)
        pure_transform(sh1, sh2, sh3, sh4, 1);

    // Results are in source_
    return size;
}


} // close namespace panache

