/*! \file
 * \brief Base class for two-body AO integrals (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/SphericalTransform.h"
#include "panache/Molecule.h"
#include "panache/TwoBodyAOInt.h"
#include "panache/AOShellCombinationsIterator.h"
#include "panache/BasisFunctionMacros.h"
#include "panache/BasisSet.h"
#include "panache/Exception.h"

namespace panache
{

static void transform2e_1(int, const SphericalTransform&, double*, double*, int);
static void transform2e_2(int, const SphericalTransform&, double*, double*, int, int, int);
static void transform2e_3(int, const SphericalTransform&, double*, double*, int, int, int);
static void transform2e_4(int, const SphericalTransform&, double*, double*, int, int);


TwoBodyAOInt::TwoBodyAOInt(const SharedBasisSet original_bs1,
                           const SharedBasisSet original_bs2,
                           const SharedBasisSet original_bs3,
                           const SharedBasisSet original_bs4)
    : original_bs1_(original_bs1),
      original_bs2_(original_bs2),
      original_bs3_(original_bs3),
      original_bs4_(original_bs4),
      curshells({{-1, -1, -1, -1}}),
      target_(0)
{
    // The derived classes allocate this memory.
    force_cartesian_ = false;
    tformbuf_ = 0;
    source_ = 0;
    natom_ = original_bs1_->molecule()->natom();  // This assumes the 4 bases come from the same molecule.
}

TwoBodyAOInt::~TwoBodyAOInt()
{
}


size_t TwoBodyAOInt::compute_shell(const AOShellCombinationsIterator& shellIter)
{
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}


double TwoBodyAOInt::compute_basisfunction(int bf1, int bf2, int bf3, int bf4)
{

    // \todo permute here so compute_shell doesn't have to?

    // first get the shell numbers
    int sh1 = original_bs1_->function_to_shell(bf1);
    int sh2 = original_bs2_->function_to_shell(bf2);
    int sh3 = original_bs3_->function_to_shell(bf3);
    int sh4 = original_bs4_->function_to_shell(bf4);

    std::array<int, 4> theseshells({{sh1, sh2, sh3, sh4}});
    if(curshells != theseshells)
    {
        // compute over the shells
        curnint = compute_shell(sh1, sh2, sh3, sh4);

        curshells = theseshells;
    }

    if(curnint == 0)
        return 0.0;

    // get only the integral we want
    const GaussianShell & gsh1 = original_bs1_->shell(sh1);
    const GaussianShell & gsh2 = original_bs2_->shell(sh2);
    const GaussianShell & gsh3 = original_bs3_->shell(sh3);
    const GaussianShell & gsh4 = original_bs4_->shell(sh4);

    //int nsh1 = gsh1.nfunction();
    int nsh2 = gsh2.nfunction();
    int nsh3 = gsh3.nfunction();
    int nsh4 = gsh4.nfunction();

    int start1 = gsh1.function_index();
    int start2 = gsh2.function_index();
    int start3 = gsh3.function_index();
    int start4 = gsh4.function_index();

    size_t index = (bf1-start1)*nsh2*nsh3*nsh4 
                 + (bf2-start2)*nsh3*nsh4
                 + (bf3-start3)*nsh4
                 + (bf4-start4);

    return buffer()[index];
}


void TwoBodyAOInt::permute_target(double *s, double *t, int sh1, int sh2, int sh3, int sh4, bool p12, bool p34, bool p13p24)
{
    const GaussianShell& s1 = bs1_->shell(sh1);
    const GaussianShell& s2 = bs2_->shell(sh2);
    const GaussianShell& s3 = bs3_->shell(sh3);
    const GaussianShell& s4 = bs4_->shell(sh4);

    int nbf1, nbf2, nbf3, nbf4;
    if(force_cartesian_)
    {
        nbf1 = s1.ncartesian();
        nbf2 = s2.ncartesian();
        nbf3 = s3.ncartesian();
        nbf4 = s4.ncartesian();
    }
    else
    {
        nbf1 = s1.nfunction();
        nbf2 = s2.nfunction();
        nbf3 = s3.nfunction();
        nbf4 = s4.nfunction();
    }

    if (!p13p24)
    {
        if (p12)
        {
            if (p34)
            {
                permute_1234_to_2143(s, t, nbf1, nbf2, nbf3, nbf4);
            }
            else
            {
                permute_1234_to_2134(s, t, nbf1, nbf2, nbf3, nbf4);
            }
        }
        else
        {
            permute_1234_to_1243(s, t, nbf1, nbf2, nbf3, nbf4);
        }
    }
    else
    {
        if (p12)
        {
            if (p34)
            {
                permute_1234_to_4321(s, t, nbf1, nbf2, nbf3, nbf4);
            }
            else
            {
                permute_1234_to_4312(s, t, nbf1, nbf2, nbf3, nbf4);
            }
        }
        else
        {
            if (p34)
            {
                permute_1234_to_3421(s, t, nbf1, nbf2, nbf3, nbf4);
            }
            else
            {
                permute_1234_to_3412(s, t, nbf1, nbf2, nbf3, nbf4);
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_1243(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf1;
    int f2=nbf2;
    int f3=nbf4;
    int f4=nbf3;
    for (int bf1=0; bf1<f1; bf1++)
    {
        for (int bf2=0; bf2<f2; bf2++)
        {
            for (int bf4=0; bf4<f4; bf4++)
            {
                for (int bf3=0; bf3<f3; bf3++)
                {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_2134(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf2;
    int f2=nbf1;
    int f3=nbf3;
    int f4=nbf4;
    for (int bf2=0; bf2<f2; bf2++)
    {
        for (int bf1=0; bf1<f1; bf1++)
        {
            for (int bf3=0; bf3<f3; bf3++)
            {
                for (int bf4=0; bf4<f4; bf4++)
                {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_2143(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf2;
    int f2=nbf1;
    int f3=nbf4;
    int f4=nbf3;
    for (int bf2=0; bf2<f2; bf2++)
    {
        for (int bf1=0; bf1<f1; bf1++)
        {
            for (int bf4=0; bf4<f4; bf4++)
            {
                for (int bf3=0; bf3<f3; bf3++)
                {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_3412(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf3;
    int f2=nbf4;
    int f3=nbf1;
    int f4=nbf2;
    for (int bf3=0; bf3<f3; bf3++)
    {
        for (int bf4=0; bf4<f4; bf4++)
        {
            for (int bf1=0; bf1<f1; bf1++)
            {
                for (int bf2=0; bf2<f2; bf2++)
                {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_4312(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf4;
    int f2=nbf3;
    int f3=nbf1;
    int f4=nbf2;
    for (int bf3=0; bf3<f3; bf3++)
    {
        for (int bf4=0; bf4<f4; bf4++)
        {
            for (int bf2=0; bf2<f2; bf2++)
            {
                for (int bf1=0; bf1<f1; bf1++)
                {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_3421(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf3;
    int f2=nbf4;
    int f3=nbf2;
    int f4=nbf1;
    for (int bf4=0; bf4<f4; bf4++)
    {
        for (int bf3=0; bf3<f3; bf3++)
        {
            for (int bf1=0; bf1<f1; bf1++)
            {
                for (int bf2=0; bf2<f2; bf2++)
                {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_4321(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4)
{
    int f1=nbf4;
    int f2=nbf3;
    int f3=nbf2;
    int f4=nbf1;
    for (int bf4=0; bf4<f4; bf4++)
    {
        for (int bf3=0; bf3<f3; bf3++)
        {
            for (int bf2=0; bf2<f2; bf2++)
            {
                for (int bf1=0; bf1<f1; bf1++)
                {
                    double *t_ptr = t + ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::pure_transform(int sh1, int sh2, int sh3, int sh4, int nchunk)
{
    const GaussianShell& s1 = bs1_->shell(sh1);
    const GaussianShell& s2 = bs2_->shell(sh2);
    const GaussianShell& s3 = bs3_->shell(sh3);
    const GaussianShell& s4 = bs4_->shell(sh4);

    // Get the transforms from the basis set
    SphericalTransform trans1(SphericalTransform::Generate(s1.am()));
    SphericalTransform trans2(SphericalTransform::Generate(s2.am()));
    SphericalTransform trans3(SphericalTransform::Generate(s3.am()));
    SphericalTransform trans4(SphericalTransform::Generate(s4.am()));

    // Get the angular momentum for each shell
    int am1 = s1.am();
    int am2 = s2.am();
    int am3 = s3.am();
    int am4 = s4.am();

    // Get number of Cartesian functions for each shell
    int nao1 = s1.ncartesian();
    int nao2 = s2.ncartesian();
    int nao3 = s3.ncartesian();
    int nao4 = s4.ncartesian();

    int nbf1 = s1.nfunction();
    int nbf2 = s2.nfunction();
    int nbf3 = s3.nfunction();
    int nbf4 = s4.nfunction();

    // Get if each shell has pure functions
    bool is_pure1 = s1.is_pure();
    bool is_pure2 = s2.is_pure();
    bool is_pure3 = s3.is_pure();
    bool is_pure4 = s4.is_pure();

    for (int ichunk=0; ichunk < nchunk; ++ichunk)
    {
        // Compute the offset in source_, and target
        size_t sourcechunkoffset = ichunk * (nao1 * nao2 * nao3 * nao4);

        // assignments to make compilers happy
        double *source1 = nullptr, *target1 = nullptr;
        double *source2 = nullptr, *target2 = nullptr;
        double *source3 = nullptr, *target3 = nullptr;
        double *source4 = nullptr, *target4 = nullptr;
        double *source = source_+sourcechunkoffset;
        double *target = target_+sourcechunkoffset;
        double *tmpbuf = tformbuf_;

        int transform_index = 8*is_pure1 + 4*is_pure2 + 2*is_pure3 + is_pure4;
        switch (transform_index)
        {
        case 0:
            break;

        case 1:
            source4 = source;
            target4 = target;
            break;

        case 2:
            source3 = source;
            target3 = target;
            break;

        case 3:
            source4 = source;
            target4 = tmpbuf;
            source3 = tmpbuf;
            target3 = target;
            break;

        case 4:
            source2 = source;
            target2 = target;
            break;

        case 5:
            source4 = source;
            target4 = tmpbuf;
            source2 = tmpbuf;
            target2 = target;
            break;

        case 6:
            source3 = source;
            target3 = tmpbuf;
            source2 = tmpbuf;
            target2 = target;
            break;

        case 7:
            source4 = source;
            target4 = tmpbuf;
            source3 = tmpbuf;
            target3 = source;
            source2 = source;
            target2 = target;
            break;

        case 8:
            source1 = source;
            target1 = target;
            break;

        case 9:
            source4 = source;
            target4 = tmpbuf;
            source1 = tmpbuf;
            target1 = target;
            break;

        case 10:
            source3 = source;
            target3 = tmpbuf;
            source1 = tmpbuf;
            target1 = target;
            break;

        case 11:
            source4 = source;
            target4 = tmpbuf;
            source3 = tmpbuf;
            target3 = source;
            source1 = source;
            target1 = target;
            break;

        case 12:
            source2 = source;
            target2 = tmpbuf;
            source1 = tmpbuf;
            target1 = target;
            break;

        case 13:
            source4 = source;
            target4 = tmpbuf;
            source2 = tmpbuf;
            target2 = source;
            source1 = source;
            target1 = target;
            break;

        case 14:
            source3 = source;
            target3 = tmpbuf;
            source2 = tmpbuf;
            target2 = source;
            source1 = source;
            target1 = target;
            break;

        case 15:
            source4 = source;
            target4 = tmpbuf;
            source3 = tmpbuf;
            target3 = source;
            source2 = source;
            target2 = tmpbuf;
            source1 = tmpbuf;
            target1 = target;
            break;
        default:
            throw RuntimeError("Invalid transform_index");  // mostly to make compilers happy
        }

        size_t size = nbf1 * nbf2 * nbf3 * nbf4;
        if (is_pure4)
        {
            transform2e_4(am4, trans4, source4, target4, nao1*nao2*nao3,nao4);
        }
        if (is_pure3)
        {
            transform2e_3(am3, trans3, source3, target3, nao1*nao2,nao3,nbf4);
        }
        if (is_pure2)
        {
            transform2e_2(am2, trans2, source2, target2, nao1,nao2,nbf3*nbf4);
        }
        if (is_pure1)
        {
            transform2e_1(am1, trans1, source1, target1, nbf2*nbf3*nbf4);
        }

        // The permute indices routines depend on the integrals being in source_
        if (is_pure1 || is_pure2 || is_pure3 || is_pure4)
            std::copy(target, target + size, source);
    }
}


/////////////////////
// HELPER FUNCTIONS
/////////////////////
static void transform2e_1(int am, const SphericalTransform& sti, double *s, double *t, int njkl)
{
    std::fill(t, t + INT_NPURE(am)*njkl, 0);

    for (auto it = sti.cbegin(); it != sti.cend(); ++it)
    {
        double *sptr = s + it->cartindex*njkl;
        double *tptr = t + it->pureindex*njkl;
        double coef = it->coef;

//        output::printf("2e_1: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        for(int jkl=0; jkl<njkl; jkl++)
            *(tptr++) += coef **(sptr++);
    }
}

static void transform2e_2(int am, const SphericalTransform& sti, double *s, double *t, int ni, int nj, int nkl)
{
    int sj = INT_NPURE(am);
    const int sjkl = nj*nkl;
    const int tjkl = sj*nkl;

    std::fill(t, t+ni*tjkl, 0);

    for (auto it = sti.cbegin(); it != sti.cend(); ++it)
    {
        double *sptr = s + it->cartindex*nkl;
        double *tptr = t + it->pureindex*nkl;
        double coef = it->coef;

//        output::printf("2e_2: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        for(int i=0; i<ni; i++,sptr+=sjkl,tptr+=tjkl)
        {
            for(int kl=0; kl<nkl; kl++)
                tptr[kl] += coef * sptr[kl];
        }
    }
}

static void transform2e_3(int am, const SphericalTransform& sti, double *s, double *t, int nij, int nk, int nl)
{
    int sk = INT_NPURE(am);
    const int skl = nk*nl;
    const int tkl = sk*nl;

    std::fill(t, t+nij*tkl, 0);

    for (auto it = sti.cbegin(); it != sti.cend(); ++it)
    {
        double *sptr = s + it->cartindex*nl;
        double *tptr = t + it->pureindex*nl;
        //output::printf("cartindex = %d, pureindex = %d\n", sti.cartindex(), sti.pureindex());

//        output::printf("2e_3: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        double coef = it->coef;
        for(int ij=0; ij<nij; ij++,sptr+=skl,tptr+=tkl)
        {
            for(int l=0; l<nl; l++)
                tptr[l] += coef * sptr[l];
        }
    }
}

// am => angular momentum of l
// sti => spherical tranformation iterator
// s => source integrals buffer
// t => target buffer
// nijk => how many i, j, k combinations are there?
// nl => how man l's are there?
static void transform2e_4(int am, const SphericalTransform& sti, double *s, double *t, int nijk, int nl)
{
    // Protect ourselves
    const int sl = nl;
    const int tl = INT_NPURE(am);

    // Clear out target memory
    std::fill(t, t+nijk*tl, 0);

    for (auto it = sti.cbegin(); it != sti.cend(); ++it)
    {
        // Starting point in source and target buffers
        double *sptr = s + it->cartindex;
        double *tptr = t + it->pureindex;

//        output::printf("2e_4: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(), sti.coef());

        // What's the coefficient we're using
        double coef = it->coef;
        for (int ijk=0; ijk<nijk; ++ijk)
        {
            // Add contribution of the source to the target
            *(tptr) += coef **(sptr);

            // skip ahead to the next ijk
            sptr += sl;
            tptr += tl;
        }
    }
}



} // close namespace panache

