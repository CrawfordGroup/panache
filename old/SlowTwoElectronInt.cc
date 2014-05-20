#include <cstring> // memset, etc
#include "SlowTwoElectronInt.h"
#include "BasisSet.h"
#include "BasisFunctionMacros.h"
#include "AOShellCombinationsIterator.h"
#include "libciomr.h"
#include "CartesianIter.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) > (b) ? (b) : (a))

namespace panache
{

namespace
{


static inline int ioff(int i)
{
    return (i*(i+1))/2;
}

} // close anonymous namespace

SlowTwoElectronInt::SlowTwoElectronInt(const std::shared_ptr<BasisSet> bs1,
                                       const std::shared_ptr<BasisSet> bs2,
                                       const std::shared_ptr<BasisSet> bs3,
                                       const std::shared_ptr<BasisSet> bs4)
    : TwoBodyAOInt(bs1,bs2,bs3,bs4)
{
    // Figure out some information to initialize libint with
    // 1. Maximum angular momentum
    int max_am = MAX(MAX(basis1()->max_am(), basis2()->max_am()), MAX(basis3()->max_am(), basis4()->max_am()));
    // 2. Maximum number of primitive combinations
    int max_nprim = basis1()->max_nprimitive() * basis2()->max_nprimitive() * basis3()->max_nprimitive() * basis4()->max_nprimitive();
    // 3. Maximum Cartesian class size
    max_cart_ = ioff(basis1()->max_am()+1) * ioff(basis2()->max_am()+1) * ioff(basis3()->max_am()+1) * ioff(basis4()->max_am()+1);

    size_t size = INT_NCART(basis1()->max_am()) * INT_NCART(basis2()->max_am()) *
                  INT_NCART(basis3()->max_am()) * INT_NCART(basis4()->max_am());

    // Used in pure_transform
    tformbuf_ = new double[size];
    memset(tformbuf_, 0, sizeof(double)*size);

    target_ = new double[size];
    memset(target_, 0, sizeof(double)*size);

    source_ = new double[size];
    memset(source_, 0, sizeof(double)*size);

    df = init_array(2 * 100);
    df[0] = 1.0;
    df[1] = 1.0;
    df[2] = 1.0;
    for (int i = 3; i < 100 * 2; i++)
    {
        df[i] = (i - 1) * df[i - 2];
    }

    fac = init_array(100);
    fac[0] = 1.0;
    for (int i = 1; i < 100; i++)
        fac[i] = fac[i - 1] * i;
    bc = block_matrix(100, 100);
    for (int i = 0; i < 100; i++)
        for (int j = 0; j <= i; j++)
            bc[i][j] = fac[i] / (fac[i - j] * fac[j]);


}

SlowTwoElectronInt::~SlowTwoElectronInt()
{
    delete[] tformbuf_;
    delete[] target_;
    delete[] source_;
    delete [] df;
    delete [] fac;
    free_block(bc);
}

size_t SlowTwoElectronInt::compute_shell(const AOShellCombinationsIterator& shellIter)
{
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}

size_t SlowTwoElectronInt::compute_shell(int sh1, int sh2, int sh3, int sh4)
{
#ifdef MINTS_TIMER
    timer_on("ERI::compute_shell");
#endif
    // Need to ensure the ordering asked by the user is valid for libint
    // compute_quartet does NOT check this. SEGFAULTS should occur if order
    // is not guaranteed.
#ifdef MINTS_TIMER
    timer_on("reorder");
#endif

    int s1, s2, s3, s4, c1, c2, c3, c4;
    int am1, am2, am3, am4, temp;
    shared_ptr<BasisSet> bs_temp;

    // AM used for ordering
    am1 = original_bs1_->shell(sh1).am();
    am2 = original_bs2_->shell(sh2).am();
    am3 = original_bs3_->shell(sh3).am();
    am4 = original_bs4_->shell(sh4).am();
    temp = am1+am2+am3+am4;

    c1 = original_bs1_->shell(sh1).ncenter();
    c2 = original_bs1_->shell(sh2).ncenter();
    c3 = original_bs1_->shell(sh3).ncenter();
    c4 = original_bs1_->shell(sh4).ncenter();

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

    size_t ncomputed = compute_quartet(sh1, sh2, sh3, sh4);

    if (ncomputed)
    {
        // Only do the following if we did any work.
        // copy the integrals to the target_
        memcpy(target_, source_, n1 * n2 * n3 * n4 *sizeof(double));
    }

#ifdef MINTS_TIMER
    timer_off("ERI::compute_shell");
#endif
    return ncomputed;
}

size_t SlowTwoElectronInt::compute_quartet(int sh1, int sh2, int sh3, int sh4)
{
#ifdef MINTS_TIMER
    timer_on("setup");
#endif

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


    CartesianIter cit1(am1);
    CartesianIter cit2(am2);
    CartesianIter cit3(am3);
    CartesianIter cit4(am4);
    cit1.start();
    cit2.start();
    cit3.start();
    cit4.start();

    int ncart1 = INT_NCART(am1);
    int ncart2 = INT_NCART(am2);
    int ncart3 = INT_NCART(am3);
    int ncart4 = INT_NCART(am4);
    size_t size = ncart1 * ncart2 * ncart3 * ncart4;


    std::array<int, 3> exp1, exp2, exp3, exp4;

    int ijkl = 0;

    for(int i = 0; i < ncart1; i++)
    {
        exp1 = cit1.exponents();

        for(int j = 0; j < ncart2; j++)
        {
            exp2 = cit2.exponents();
            for(int k = 0; k < ncart3; k++)
            {
                exp3 = cit3.exponents();
                for(int l = 0; l < ncart4; l++)
                {
                    exp4 = cit4.exponents();

                    source_[ijkl] = 0;

                    // Loop over primitives
                    for(int a = 0; a < s1.nprimitive(); a++)
                        for(int b = 0; b < s2.nprimitive(); b++)
                            for(int c = 0; c < s3.nprimitive(); c++)
                                for(int d = 0; d < s4.nprimitive(); d++)
                                    source_[ijkl] += eri(
                                                         exp1[0], exp1[1], exp1[2], s1.exp(a), A,
                                                         exp2[0], exp2[1], exp2[2], s2.exp(a), B,
                                                         exp3[0], exp3[1], exp3[2], s3.exp(a), C,
                                                         exp4[0], exp4[1], exp4[2], s4.exp(a), D,
                                                         0)*s1.coef(a)*s2.coef(b)*s3.coef(c)*s4.coef(d);
                }

            }

        }

    }

    // Transform the integrals into pure angular momentum
    if (!force_cartesian_)
        pure_transform(sh1, sh2, sh3, sh4, 1);

    // Results are in source_
    return size;
}



double SlowTwoElectronInt::eri(unsigned int l1, unsigned int m1, unsigned int n1, double alpha1,
                               const double* A, unsigned int l2, unsigned int m2, unsigned int n2,
                               double alpha2, const double* B, unsigned int l3, unsigned int m3,
                               unsigned int n3, double alpha3, const double* C, unsigned int l4,
                               unsigned int m4, unsigned int n4, double alpha4, const double* D,
                               int norm_flag)
{

    const double gammap = alpha1 + alpha2;
    const double Px = (alpha1 * A[0] + alpha2 * B[0]) / gammap;
    const double Py = (alpha1 * A[1] + alpha2 * B[1]) / gammap;
    const double Pz = (alpha1 * A[2] + alpha2 * B[2]) / gammap;
    const double PAx = Px - A[0];
    const double PAy = Py - A[1];
    const double PAz = Pz - A[2];
    const double PBx = Px - B[0];
    const double PBy = Py - B[1];
    const double PBz = Pz - B[2];
    const double AB2 = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1]
                       - B[1]) + (A[2] - B[2]) * (A[2] - B[2]);

    const double gammaq = alpha3 + alpha4;
    const double gammapq = gammap * gammaq / (gammap + gammaq);
    const double Qx = (alpha3 * C[0] + alpha4 * D[0]) / gammaq;
    const double Qy = (alpha3 * C[1] + alpha4 * D[1]) / gammaq;
    const double Qz = (alpha3 * C[2] + alpha4 * D[2]) / gammaq;
    const double QCx = Qx - C[0];
    const double QCy = Qy - C[1];
    const double QCz = Qz - C[2];
    const double QDx = Qx - D[0];
    const double QDy = Qy - D[1];
    const double QDz = Qz - D[2];
    const double CD2 = (C[0] - D[0]) * (C[0] - D[0]) + (C[1] - D[1]) * (C[1]
                       - D[1]) + (C[2] - D[2]) * (C[2] - D[2]);

    const double PQx = Px - Qx;
    const double PQy = Py - Qy;
    const double PQz = Pz - Qz;
    const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

    int u1, u2, v1, v2, w1, w2, tx, ty, tz, txmax, tymax, tzmax;
    int i, j, k;
    int lp, lq, mp, mq, np, nq;
    int zeta;
    double *flp, *flq, *fmp, *fmq, *fnp, *fnq;
    double *F;
    double K1, K2;
    double Gx, Gy, Gz;
    double pfac;
    double result = 0.0;
    double tmp;
    int u1max, u2max, v1max, v2max, w1max, w2max;

    K1 = exp(-alpha1 * alpha2 * AB2 / gammap);
    K2 = exp(-alpha3 * alpha4 * CD2 / gammaq);
    pfac = 2 * std::pow(M_PI, 2.5) * K1 * K2 / (gammap * gammaq
            * sqrt(gammap + gammaq));

    if (fac == NULL)
    {
        fac = init_array(100);
        fac[0] = 1.0;
        for (i = 1; i < 100; i++)
            fac[i] = fac[i - 1] * i;
        bc = block_matrix(100, 100);
        for (i = 0; i < 100; i++)
            for (j = 0; j <= i; j++)
                bc[i][j] = fac[i] / (fac[i - j] * fac[j]);
    }

    if (norm_flag > 0)
    {
        pfac *= norm_const(l1, m1, n1, alpha1, A);
        pfac *= norm_const(l2, m2, n2, alpha2, B);
        pfac *= norm_const(l3, m3, n3, alpha3, C);
        pfac *= norm_const(l4, m4, n4, alpha4, D);
    }

    F = init_array(l1 + l2 + l3 + l4 + m1 + m2 + m3 + m4 + n1 + n2 + n3 + n4 + 1);
    calc_f(F, l1 + l2 + l3 + l4 + m1 + m2 + m3 + m4 + n1 + n2 + n3 + n4,
           PQ2 * gammapq);

    flp = init_array(l1 + l2 + 1);
    for (k = 0; k <= l1 + l2; k++)
        for (i = 0; i <= MIN(k,l1); i++)
        {
            j = k - i;
            if (j > l2)
                continue;
            tmp = bc[l1][i] * bc[l2][j];
            if (l1 - i > 0)
                tmp *= pow(PAx, l1 - i);
            if (l2 - j > 0)
                tmp *= pow(PBx, l2 - j);
            flp[k] += tmp;
        }
    fmp = init_array(m1 + m2 + 1);
    for (k = 0; k <= m1 + m2; k++)
        for (i = 0; i <= MIN(k,m1); i++)
        {
            j = k - i;
            if (j > m2)
                continue;
            tmp = bc[m1][i] * bc[m2][j];
            if (m1 - i > 0)
                tmp *= pow(PAy, m1 - i);
            if (m2 - j > 0)
                tmp *= pow(PBy, m2 - j);
            fmp[k] += tmp;
        }
    fnp = init_array(n1 + n2 + 1);
    for (k = 0; k <= n1 + n2; k++)
        for (i = 0; i <= MIN(k,n1); i++)
        {
            j = k - i;
            if (j > n2)
                continue;
            tmp = bc[n1][i] * bc[n2][j];
            if (n1 - i > 0)
                tmp *= pow(PAz, n1 - i);
            if (n2 - j > 0)
                tmp *= pow(PBz, n2 - j);
            fnp[k] += tmp;
        }
    flq = init_array(l3 + l4 + 1);
    for (k = 0; k <= l3 + l4; k++)
        for (i = 0; i <= MIN(k,l3); i++)
        {
            j = k - i;
            if (j > l4)
                continue;
            tmp = bc[l3][i] * bc[l4][j];
            if (l3 - i > 0)
                tmp *= pow(QCx, l3 - i);
            if (l4 - j > 0)
                tmp *= pow(QDx, l4 - j);
            flq[k] += tmp;
        }
    fmq = init_array(m3 + m4 + 1);
    for (k = 0; k <= m3 + m4; k++)
        for (i = 0; i <= MIN(k,m3); i++)
        {
            j = k - i;
            if (j > m4)
                continue;
            tmp = bc[m3][i] * bc[m4][j];
            if (m3 - i > 0)
                tmp *= pow(QCy, m3 - i);
            if (m4 - j > 0)
                tmp *= pow(QDy, m4 - j);
            fmq[k] += tmp;
        }
    fnq = init_array(n3 + n4 + 1);
    for (k = 0; k <= n3 + n4; k++)
        for (i = 0; i <= MIN(k,n3); i++)
        {
            j = k - i;
            if (j > n4)
                continue;
            tmp = bc[n3][i] * bc[n4][j];
            if (n3 - i > 0)
                tmp *= pow(QCz, n3 - i);
            if (n4 - j > 0)
                tmp *= pow(QDz, n4 - j);
            fnq[k] += tmp;
        }

    for (lp = 0; lp <= l1 + l2; lp++)
        for (lq = 0; lq <= l3 + l4; lq++)
        {
            u1max = lp / 2;
            u2max = lq / 2;
            for (u1 = 0; u1 <= u1max; u1++)
                for (u2 = 0; u2 <= u2max; u2++)
                {
                    Gx = pow(-1, lp) * flp[lp] * flq[lq] * fac[lp] * fac[lq]
                         * pow(gammap, u1 - lp) * pow(gammaq, u2 - lq) * fac[lp + lq - 2
                                 * u1 - 2 * u2] * pow(gammapq, lp + lq - 2 * u1 - 2 * u2)
                         / (fac[u1] * fac[u2] * fac[lp - 2 * u1] * fac[lq - 2 * u2]);
                    for (mp = 0; mp <= m1 + m2; mp++)
                        for (mq = 0; mq <= m3 + m4; mq++)
                        {
                            v1max = mp / 2;
                            v2max = mq / 2;
                            for (v1 = 0; v1 <= v1max; v1++)
                                for (v2 = 0; v2 <= v2max; v2++)
                                {
                                    Gy = pow(-1, mp) * fmp[mp] * fmq[mq] * fac[mp] * fac[mq]
                                         * pow(gammap, v1 - mp) * pow(gammaq, v2 - mq) * fac[mp
                                                 + mq - 2 * v1 - 2 * v2] * pow(gammapq,
                                                         mp + mq - 2 * v1 - 2 * v2)
                                         / (fac[v1] * fac[v2] * fac[mp - 2 * v1]
                                            * fac[mq - 2 * v2]);
                                    for (np = 0; np <= n1 + n2; np++)
                                        for (nq = 0; nq <= n3 + n4; nq++)
                                        {
                                            w1max = np / 2;
                                            w2max = nq / 2;
                                            for (w1 = 0; w1 <= w1max; w1++)
                                                for (w2 = 0; w2 <= w2max; w2++)
                                                {
                                                    Gz = pow(-1, np) * fnp[np] * fnq[nq] * fac[np]
                                                         * fac[nq] * pow(gammap, w1 - np) * pow(gammaq,
                                                                 w2 - nq)
                                                         * fac[np + nq - 2 * w1 - 2 * w2]
                                                         * pow(gammapq, np + nq - 2 * w1 - 2 * w2)
                                                         / (fac[w1] * fac[w2] * fac[np - 2 * w1] * fac[nq
                                                                 - 2 * w2]);
                                                    txmax = (lp + lq - 2 * u1 - 2 * u2) / 2;
                                                    tymax = (mp + mq - 2 * v1 - 2 * v2) / 2;
                                                    tzmax = (np + nq - 2 * w1 - 2 * w2) / 2;
                                                    for (tx = 0; tx <= txmax; tx++)
                                                        for (ty = 0; ty <= tymax; ty++)
                                                            for (tz = 0; tz <= tzmax; tz++)
                                                            {
                                                                zeta = lp + lq + mp + mq + np + nq - 2 * u1 - 2
                                                                       * u2 - 2 * v1 - 2 * v2 - 2 * w1 - 2 * w2
                                                                       - tx - ty - tz;
                                                                result += Gx * Gy * Gz * F[zeta]
                                                                          * pow(-1, tx + ty + tz) * pow(
                                                                              PQx,
                                                                              lp + lq - 2
                                                                              * u1 - 2
                                                                              * u2 - 2
                                                                              * tx)
                                                                          * pow(PQy,
                                                                                mp + mq - 2 * v1 - 2 * v2 - 2 * ty)
                                                                          * pow(PQz,
                                                                                np + nq - 2 * w1 - 2 * w2 - 2 * tz)
                                                                          / (pow(
                                                                                 4,
                                                                                 u1 + u2 + tx + v1 + v2 + ty + w1
                                                                                 + w2 + tz) * pow(gammapq, tx)
                                                                             * pow(gammapq, ty) * pow(gammapq, tz)
                                                                             * fac[lp + lq - 2 * u1 - 2 * u2 - 2
                                                                                   * tx] * fac[tx] * fac[mp + mq - 2
                                                                                           * v1 - 2 * v2 - 2 * ty] * fac[ty]
                                                                             * fac[np + nq - 2 * w1 - 2 * w2 - 2
                                                                                   * tz] * fac[tz]);
                                                            }
                                                }
                                        }
                                }
                        }
                }
        }

    free_array(F);
    free_array(flp);
    free_array(fmp);
    free_array(fnp);
    free_array(flq);
    free_array(fmq);
    free_array(fnq);

    return result * pfac;
}


/*!
 calc_f()

 This function computes infamous integral Fn(t). For its definition
 see Obara and Saika paper, or Shavitt's chapter in the
 Methods in Computational Physics book (see reference below).
 This piece of code is from Dr. Justin Fermann's program CINTS

 \ingroup (QT)
 */
void SlowTwoElectronInt::calc_f(double *F, int n, double t)
{
    int i, m, k;
    int m2;
    double t2;
    double num;
    double sum;
    double term1, term2;
    static double K = 1.0 / M_2_SQRTPI;
    double et;

    if (df == NULL)
    {
        df = init_array(2 * 100);
        df[0] = 1.0;
        df[1] = 1.0;
        df[2] = 1.0;
        for (i = 3; i < 100 * 2; i++)
        {
            df[i] = (i - 1) * df[i - 2];
        }
    }

    if (t > 20.0)   /* For big t's do upward recursion */
    {
        t2 = 2 * t;
        et = exp(-t);
        t = sqrt(t);
        F[0] = K * erf(t) / t;
        for (m = 0; m <= n - 1; m++)
        {
            F[m + 1] = ((2 * m + 1) * F[m] - et) / (t2);
        }
    }
    else
    {
        /* For smaller t's compute F with highest n using
         asymptotic series (see I. Shavitt in
         Methods in Computational Physics, ed. B. Alder eta l,
         vol 2, 1963, page 8) */
        et = exp(-t);
        t2 = 2 * t;
        m2 = 2 * n;
        num = df[m2];
        i = 0;
        sum = 1.0 / (m2 + 1);
        do
        {
            i++;
            num = num * t2;
            term1 = num / df[m2 + 2 * i + 2];
            sum += term1;
        }
        while (fabs(term1) > 1.0E-17 && i < 100);
        F[n] = sum * et;
        for (m = n - 1; m >= 0; m--)   /* And then do downward recursion */
        {
            F[m] = (t2 * F[m + 1] + et) / (2 * m + 1);
        }
    }
}

/*!
 norm_const()

 \ingroup (QT)
 */
double SlowTwoElectronInt::norm_const(unsigned int l1, unsigned int m1, unsigned int n1,
                                      double alpha1, const double* A)
{
    return pow(2 * alpha1 / M_PI, 0.75) * pow(4 * alpha1, 0.5 * (l1 + m1 + n1))
           / sqrt(df[2 * l1] * df[2 * m1] * df[2 * n1]);
}


} // close namespace panache


