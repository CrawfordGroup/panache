#include <cmath>
#include "LibintTwoElectronInt.h"
#include "BasisSet.h"
#include "BasisFunctionMacros.h"
#include "AOShellCombinationsIterator.h"
#include "Fjt.h"
#include "Lapack.h"
#include "PhysConst.h"
#include "Exception.h"

namespace panache
{

namespace
{

static inline int ioff(int i)
{
    return (i*(i+1))/2;
}

} // end anonymous namespace



LibintTwoElectronInt::LibintTwoElectronInt(const SharedBasisSet bs1,
                                           const SharedBasisSet bs2,
                                           const SharedBasisSet bs3,
                                           const SharedBasisSet bs4)
    : TwoBodyAOInt(bs1,bs2,bs3,bs4)
{
    // Initialize libint static data
    init_libint_base();

    // Figure out some information to initialize libint with
    // 1. Maximum angular momentum
    int max_am = std::max(std::max(basis1()->max_am(), basis2()->max_am()), std::max(basis3()->max_am(), basis4()->max_am()));
    // 2. Maximum number of primitive combinations
    int max_nprim = basis1()->max_nprimitive() * basis2()->max_nprimitive() * basis3()->max_nprimitive() * basis4()->max_nprimitive();
    // 3. Maximum Cartesian class size
    max_cart_ = ioff(basis1()->max_am()+1) * ioff(basis2()->max_am()+1) * ioff(basis3()->max_am()+1) * ioff(basis4()->max_am()+1);

    // Make sure libint is compiled to handle our max AM
    if (max_am >= LIBINT_MAX_AM)
    {
        throw RuntimeError("ERI - libint cannot handle angular momentum this high.\n"
                           "In a fresh object directory, reconfigure libint for higher angular momentum, then recompile.");
    }

    try
    {
        // Initialize libint
        init_libint(&libint_, max_am, max_nprim);
    }
    catch (std::bad_alloc& e)
    {
        throw RuntimeError("Error allocating memory for libint.\n");
    }
    size_t size = INT_NCART(basis1()->max_am()) * INT_NCART(basis2()->max_am()) *
                  INT_NCART(basis3()->max_am()) * INT_NCART(basis4()->max_am());

    // Used in pure_transform
    tformbuf_ = new double[size];
    std::fill(tformbuf_, tformbuf_+size, 0);

    target_ = new double[size];
    std::fill(target_, target_+size, 0);

    source_ = new double[size];
    std::fill(source_, source_+size, 0);
}

LibintTwoElectronInt::~LibintTwoElectronInt()
{
    delete[] tformbuf_;
    delete[] target_;
    delete[] source_;
    free_libint(&libint_);
}

size_t LibintTwoElectronInt::compute_shell(const AOShellCombinationsIterator& shellIter)
{
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}

size_t LibintTwoElectronInt::compute_shell(int sh1, int sh2, int sh3, int sh4)
{
    // Need to ensure the ordering asked by the user is valid for libint
    // compute_quartet does NOT check this. SEGFAULTS should occur if order
    // is not guaranteed.

    int s1, s2, s3, s4;
    int am1, am2, am3, am4, temp;
    shared_ptr<BasisSet> bs_temp;

    p13p24_ = false;
    p12_ = false;
    p34_ = false;

    // AM used for ordering
    am1 = original_bs1_->shell(sh1).am();
    am2 = original_bs2_->shell(sh2).am();
    am3 = original_bs3_->shell(sh3).am();
    am4 = original_bs4_->shell(sh4).am();
    temp = am1+am2+am3+am4;

    // TODO: Check this!
//	if (c1 == c2 && c1 == c3 && c1 && c4 && temp % 2 != 0) {
//		return 0;
//	}

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

    // Save the original requested shell ordering. The pre-computed shell pair information
    // requires the original ordering.
    osh1_ = sh1;
    osh2_ = sh2;
    osh3_ = sh3;
    osh4_ = sh4;

    // l(a) >= l(b), l(c) >= l(d), and l(c) + l(d) >= l(a) + l(b).
    if (am1 >= am2)
    {
        s1 = sh1;
        s2 = sh2;

        bs1_ = original_bs1_;
        bs2_ = original_bs2_;
    }
    else
    {
        s1 = sh2;
        s2 = sh1;

        bs1_ = original_bs2_;
        bs2_ = original_bs1_;

        p12_ = true;
    }

    if (am3 >= am4)
    {
        s3 = sh3;
        s4 = sh4;

        bs3_ = original_bs3_;
        bs4_ = original_bs4_;

    }
    else
    {
        s3 = sh4;
        s4 = sh3;

        bs3_ = original_bs4_;
        bs4_ = original_bs3_;

        p34_ = true;
    }

    if ((am1 + am2) > (am3 + am4))
    {
        // Swap s1 and s2 with s3 and s4
        temp = s1;
        s1 = s3;
        s3 = temp;

        temp = s2;
        s2 = s4;
        s4 = temp;

        bs_temp = bs1_;
        bs1_ = bs3_;
        bs3_ = bs_temp;

        bs_temp = bs2_;
        bs2_ = bs4_;
        bs4_ = bs_temp;

        p13p24_ = true;
    }

    // s1, s2, s3, s4 contain the shells to do in libint order
    size_t ncomputed = compute_quartet(s1, s2, s3, s4);
    if (ncomputed)
    {
        // Only do the following if we did any work.

        // Permute integrals back, if needed
        if (p12_ || p34_ || p13p24_)
        {
            permute_target(source_, target_, s1, s2, s3, s4, p12_, p34_, p13p24_);
        }
        else
        {
            // copy the integrals to the target_
            std::copy(source_, source_ + n1 * n2 * n3 * n4, target_);
        }
    }

    return ncomputed;
}

size_t LibintTwoElectronInt::compute_quartet(int sh1, int sh2, int sh3, int sh4)
{

    const GaussianShell& s1 = bs1_->shell(sh1);
    const GaussianShell& s2 = bs2_->shell(sh2);
    const GaussianShell& s3 = bs3_->shell(sh3);
    const GaussianShell& s4 = bs4_->shell(sh4);

    int am1 = s1.am();
    int am2 = s2.am();
    int am3 = s3.am();
    int am4 = s4.am();
    int am = am1 + am2 + am3 + am4; // total am
    int nprim1;
    int nprim2;
    int nprim3;
    int nprim4;
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

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);
    double CD2 = 0.0;
    CD2 += (C[0] - D[0]) * (C[0] - D[0]);
    CD2 += (C[1] - D[1]) * (C[1] - D[1]);
    CD2 += (C[2] - D[2]) * (C[2] - D[2]);

    libint_.AB[0] = A[0] - B[0];
    libint_.AB[1] = A[1] - B[1];
    libint_.AB[2] = A[2] - B[2];
    libint_.CD[0] = C[0] - D[0];
    libint_.CD[1] = C[1] - D[1];
    libint_.CD[2] = C[2] - D[2];



    // Prepare all the data needed by libint
    size_t nprim = 0;
    nprim1 = s1.nprimitive();
    nprim2 = s2.nprimitive();
    nprim3 = s3.nprimitive();
    nprim4 = s4.nprimitive();

    const double *a1s = s1.exps();
    const double *a2s = s2.exps();
    const double *a3s = s3.exps();
    const double *a4s = s4.exps();
    const double *c1s = s1.coefs();
    const double *c2s = s2.coefs();
    const double *c3s = s3.coefs();
    const double *c4s = s4.coefs();

    for (int p1=0; p1<nprim1; ++p1)
    {
        double a1 = a1s[p1];
        double c1 = c1s[p1];
        for (int p2=0; p2<nprim2; ++p2)
        {
            double a2 = a2s[p2];
            double c2 = c2s[p2];
            double zeta = a1 + a2;
            double ooz = 1.0/zeta;
            double oo2z = 1.0/(2.0 * zeta);

            double PA[3];
            double P[3];

            P[0] = (a1*A[0] + a2*B[0])*ooz;
            P[1] = (a1*A[1] + a2*B[1])*ooz;
            P[2] = (a1*A[2] + a2*B[2])*ooz;
            PA[0] = P[0] - A[0];
            PA[1] = P[1] - A[1];
            PA[2] = P[2] - A[2];

            double Sab = pow(M_PI*ooz, 3.0/2.0) * exp(-a1*a2*ooz*AB2) * c1 * c2;

            for (int p3=0; p3<nprim3; ++p3)
            {
                double a3 = a3s[p3];
                double c3 = c3s[p3];
                for (int p4=0; p4<nprim4; ++p4)
                {
                    double a4 = a4s[p4];
                    double c4 = c4s[p4];
                    double nu = a3 + a4;
                    double oon = 1.0/nu;
                    double oo2n = 1.0/(2.0*nu);
                    double oo2zn = 1.0/(2.0*(zeta+nu));
                    double rho = (zeta*nu)/(zeta+nu);
                    double oo2rho = 1.0 / (2.0*rho);

                    double QC[3], WP[3], WQ[3];
                    double Q[3], W[3], a3C[3], a4D[3];

                    a3C[0] = a3*C[0];
                    a3C[1] = a3*C[1];
                    a3C[2] = a3*C[2];

                    a4D[0] = a4*D[0];
                    a4D[1] = a4*D[1];
                    a4D[2] = a4*D[2];

                    Q[0] = (a3C[0] + a4D[0])*oon;
                    Q[1] = (a3C[1] + a4D[1])*oon;
                    Q[2] = (a3C[2] + a4D[2])*oon;

                    QC[0] = Q[0] - C[0];
                    QC[1] = Q[1] - C[1];
                    QC[2] = Q[2] - C[2];

                    double PQ2 = 0.0;
                    PQ2 += (P[0] - Q[0]) * (P[0] - Q[0]);
                    PQ2 += (P[1] - Q[1]) * (P[1] - Q[1]);
                    PQ2 += (P[2] - Q[2]) * (P[2] - Q[2]);

                    W[0] = (zeta*P[0] + nu*Q[0]) / (zeta + nu);
                    W[1] = (zeta*P[1] + nu*Q[1]) / (zeta + nu);
                    W[2] = (zeta*P[2] + nu*Q[2]) / (zeta + nu);
                    WP[0] = W[0] - P[0];
                    WP[1] = W[1] - P[1];
                    WP[2] = W[2] - P[2];
                    WQ[0] = W[0] - Q[0];
                    WQ[1] = W[1] - Q[1];
                    WQ[2] = W[2] - Q[2];

                    for (int i=0; i<3; ++i)
                    {
                        libint_.PrimQuartet[nprim].U[0][i] = PA[i];
                        libint_.PrimQuartet[nprim].U[2][i] = QC[i];
                        libint_.PrimQuartet[nprim].U[4][i] = WP[i];
                        libint_.PrimQuartet[nprim].U[5][i] = WQ[i];
                    }
                    libint_.PrimQuartet[nprim].oo2z = oo2z;
                    libint_.PrimQuartet[nprim].oo2n = oo2n;
                    libint_.PrimQuartet[nprim].oo2zn = oo2zn;
                    libint_.PrimQuartet[nprim].poz = rho * ooz;
                    libint_.PrimQuartet[nprim].pon = rho * oon;
                    libint_.PrimQuartet[nprim].oo2p = oo2rho;

                    double T = rho * PQ2;
                    fjt_->set_rho(rho);
                    double * restrict F = fjt_->values(am, T);

                    // Modify F to include overlap of ab and cd, eqs 14, 15, 16 of libint manual
                    double Scd = pow(M_PI*oon, 3.0/2.0) * exp(-a3*a4*oon*CD2) * c3 * c4;
                    double val = 2.0 * sqrt(rho * M_1_PI) * Sab * Scd;

                    for (int i=0; i<=am; ++i)
                    {
                        libint_.PrimQuartet[nprim].F[i] = F[i] * val;
                    }
                    nprim++;
                }
            }
        }
    }

    // How many are there?
    size_t size = INT_NCART(am1) * INT_NCART(am2) * INT_NCART(am3) * INT_NCART(am4);


    // Compute the integral
    if (am)
    {
        double *target_ints;

        target_ints = build_eri[am1][am2][am3][am4](&libint_, nprim);

        std::copy(target_ints, target_ints + size, source_);
    }
    else
    {
        // Handle (ss|ss)
        double temp = 0.0;
        for (size_t i=0; i<nprim; ++i)
            temp += (double)libint_.PrimQuartet[i].F[0];
        source_[0] = temp;
//        output::printf("s-functions = %8.5f\n", temp);
    }


    // The following two functions time themselves.

    // Normalize the integrals for angular momentum
    //normalize_am(s1, s2, s3, s4);

    // Transform the integrals into pure angular momentum
    if (!force_cartesian_)
        pure_transform(sh1, sh2, sh3, sh4, 1);

    // Results are in source_
    return size;
}

} // close namespace panache



