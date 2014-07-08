/*! \file
 * \brief Base class for libERD-based ERI (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/ERDTwoElectronInt.h"
#include "panache/BasisSet.h"
#include "panache/BasisFunctionMacros.h"
#include "panache/AOShellCombinationsIterator.h"
#include "panache/Output.h"

#define DEBUG 0

// Name mangling
#define FC_SYMBOL 2

#if FC_SYMBOL == 1
#define C_ERD__GENER_ERI_BATCH erd__gener_eri_batch
#define C_ERD__MEMORY_CSGTO erd__memory_csgto
#define C_ERD__MEMORY_ERI_BATCH erd__memory_eri_batch
#elif FC_SYMBOL == 2
#define C_ERD__GENER_ERI_BATCH erd__gener_eri_batch_
#define C_ERD__MEMORY_CSGTO erd__memory_csgto_
#define C_ERD__MEMORY_ERI_BATCH erd__memory_eri_batch_
#elif FC_SYMBOL == 3
#define C_ERD__GENER_ERI_BATCH ERD__GENER_ERI_BATCH
#define C_ERD__MEMORY_CSGTO ERD__MEMORY_CSGTO
#define C_ERD__MEMORY_ERI_BATCH ERD__MEMORY_ERI_BATCH
#elif FC_SYMBOL == 4
#define C_ERD__GENER_ERI_BATCH ERD__GENER_ERI_BATCH_
#define C_ERD__MEMORY_CSGTO ERD__MEMORY_CSGTO_
#define C_ERD__MEMORY_ERI_BATCH ERD__MEMORY_ERI_BATCH_
#else
#error FC_SYMBOL is not defined
#endif

extern "C" {
void C_ERD__GENER_ERI_BATCH(const F_INT &imax, const F_INT &zmax, const F_INT &nalpha, const F_INT &ncoeff,
                            const F_INT &ncsum, const F_INT &ncgto1, const F_INT &ncgto2,
                            const F_INT &ncgto3, const F_INT &ncgto4, const F_INT &npgto1,
                            const F_INT &npgto2, const F_INT &npgto3, const F_INT &npgto4,
                            const F_INT &shell1, const F_INT &shell2, const F_INT &shell3,
                            const F_INT &shell4, const double &x1, const double &y1, const double &z1,
                            const double &x2,const double &y2,const double &z2, const double &x3,
                            const double &y3,const double &z3,const double &x4, const double &y4, const double &z4,
                            const double *alpha, const double *cc, const F_INT *ccbeg, const F_INT *ccend,
                            const F_BOOL &spheric,  const F_BOOL &screen, F_INT *icore,
                            F_INT &nbatch, F_INT & nfirst, double *zcore );
void C_ERD__MEMORY_ERI_BATCH(const F_INT &nalpha, const F_INT &ncoeff,
                             const F_INT &ncgto1, const F_INT &ncgto2, const F_INT &ncgto3, const F_INT &ncgto4,
                             const F_INT &npgto1, const F_INT &npgto2, const F_INT &npgto3, const F_INT &npgto4,
                             const F_INT &shell1, const F_INT &shell2, const F_INT &shell3, const F_INT &shell4,
                             const double &x1, const double &y1, const double &z1, const double &x2, const double &y2,
                             const double &z2, const double &x3, const double &y3, const double &z3, const double &x4,
                             const double &y4, const double &z4, const double *alpha, const double *cc, const F_BOOL &spheric,
                             F_INT &imin, F_INT &iopt, F_INT &zmin, F_INT &zopt);
}

namespace panache {

ERDTwoElectronInt::ERDTwoElectronInt(const SharedBasisSet bs1,
               const SharedBasisSet bs2,
               const SharedBasisSet bs3,
               const SharedBasisSet bs4)
    : TwoBodyAOInt(bs1, bs2, bs3, bs4), 
      d_buffer_size_(0L),
      i_buffer_size_(0L),
      buffer_offset_(0L)
{
    bs1_ = original_bs1_;
    bs2_ = original_bs2_;
    bs3_ = original_bs3_;
    bs4_ = original_bs4_;
    same_bs_ = (bs1_ == bs2_ && bs1_ == bs3_ && bs1_ == bs4_);

    size_t max_cart = INT_NCART(basis1()->max_am()) * INT_NCART(basis2()->max_am()) *
                      INT_NCART(basis3()->max_am()) * INT_NCART(basis4()->max_am());

    tformbuf_ = new double[max_cart];
    std::fill(tformbuf_, tformbuf_+max_cart, 0);


    target_ = new double[max_cart];
    std::fill(target_, target_ + max_cart, 0);

    size_t max_nprim = basis1()->max_nprimitive() +
                       basis2()->max_nprimitive() +
                       basis3()->max_nprimitive() +
                       basis4()->max_nprimitive();
    cc_ = new double[max_nprim];
    alpha_ = new double[max_nprim];

    screen_ = 0;
    spheric_ = 0;

    // Ask ERD for the maximum amount of scratch it'll need
    compute_scratch_size();

    dscratch_ = new double[d_buffer_size_];
    std::fill(dscratch_, dscratch_+d_buffer_size_, 0);

    iscratch_ = new F_INT[i_buffer_size_];
}


ERDTwoElectronInt::~ERDTwoElectronInt()
{
    delete[] alpha_;
    delete[] cc_;
    delete[] tformbuf_;
    delete[] target_;
    delete[] dscratch_;
    delete[] iscratch_;
}


size_t ERDTwoElectronInt::compute_shell(const AOShellCombinationsIterator& shellIter)
{
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}


void ERDTwoElectronInt::compute_scratch_size()
{
    int shell1 = 0;
    int npgto1 = 0;
    for(int shell = 0; shell < original_bs1_->nshell(); ++shell)
        if(original_bs1_->shell(shell).nprimitive() > npgto1){
            npgto1 = original_bs1_->shell(shell).nprimitive();
            shell1 = shell;
        }
    int shell2 = 0;
    int npgto2 = 0;
    for(int shell = 0; shell < original_bs2_->nshell(); ++shell)
        if(original_bs2_->shell(shell).nprimitive() > npgto2){
            npgto2 = original_bs2_->shell(shell).nprimitive();
            shell2 = shell;
        }
    int shell3 = 0;
    int npgto3 = 0;
    for(int shell = 0; shell < original_bs3_->nshell(); ++shell)
        if(original_bs3_->shell(shell).nprimitive() > npgto3){
            npgto3 = original_bs3_->shell(shell).nprimitive();
            shell3 = shell;
        }
    int shell4 = 0;
    int npgto4 = 0;
    for(int shell = 0; shell < original_bs4_->nshell(); ++shell)
        if(original_bs4_->shell(shell).nprimitive() > npgto4){
            npgto4 = original_bs4_->shell(shell).nprimitive();
            shell4 = shell;
        }
    const GaussianShell &gs1 = original_bs1_->shell(shell1);
    const GaussianShell &gs2 = original_bs2_->shell(shell2);
    const GaussianShell &gs3 = original_bs3_->shell(shell3);
    const GaussianShell &gs4 = original_bs4_->shell(shell4);
    double x1 = 1.0;
    double y1 = 1.0;
    double z1 = 1.0;
    double x2 = 2.0;
    double y2 = 2.0;
    double z2 = 2.0;
    double x3 = 3.0;
    double y3 = 3.0;
    double z3 = 3.0;
    double x4 = 4.0;
    double y4 = 4.0;
    double z4 = 4.0;
    F_INT ncgto1 = 1;
    F_INT ncgto2 = 1;
    F_INT ncgto3 = 1;
    F_INT ncgto4 = 1;
    F_INT npgto = npgto1 + npgto2 + npgto3 + npgto4;
    F_INT am1 = original_bs1_->max_am();
    F_INT am2 = original_bs2_->max_am();
    F_INT am3 = original_bs3_->max_am();
    F_INT am4 = original_bs4_->max_am();
    F_INT imin = 0;
    F_INT iopt = 0;
    F_INT zmin = 0;
    F_INT zopt = 0;
    F_INT last_pgto = 0;
    for(int pgto1 = 0; pgto1 < npgto1; ++pgto1){
        cc_[last_pgto] = gs1.coef(pgto1);
        alpha_[last_pgto] = gs1.exp(pgto1);
        ++last_pgto;
    }
    for(int pgto2 = 0; pgto2 < npgto2; ++pgto2){
        cc_[last_pgto] = gs2.coef(pgto2);
        alpha_[last_pgto] = gs2.exp(pgto2);
        ++last_pgto;
    }
    for(int pgto3 = 0; pgto3 < npgto3; ++pgto3){
        cc_[last_pgto] = gs3.coef(pgto3);
        alpha_[last_pgto] = gs3.exp(pgto3);
        ++last_pgto;
    }
    for(int pgto4 = 0; pgto4 < npgto4; ++pgto4){
        cc_[last_pgto] = gs4.coef(pgto4);
        alpha_[last_pgto] = gs4.exp(pgto4);
        ++last_pgto;
    }
    // Compute the amount of memory needed for the largest quartet
    C_ERD__MEMORY_ERI_BATCH(npgto, npgto, ncgto1, ncgto2, ncgto3, ncgto4,
                            npgto1, npgto2, npgto3, npgto4, am1, am2, am3, am4,
                            x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, 
                            alpha_, cc_, spheric_, imin, iopt, zmin, zopt);
#if DEBUG
    output::printf("\timin %d iopt %d zmin %d zopt %d\n", imin, iopt, zmin, zopt);
#endif
    i_buffer_size_ = iopt;
    d_buffer_size_ = zopt;

    // We always treat basis sets as segmented, right now.
    ccbeg_[0] = 1;
    ccbeg_[1] = 1;
    ccbeg_[2] = 1;
    ccbeg_[3] = 1;
}

size_t ERDTwoElectronInt::compute_shell(int shell_i, int shell_j, int shell_k, int shell_l)
{
    const GaussianShell &gs1 = original_bs1_->shell(shell_i);
    const GaussianShell &gs2 = original_bs2_->shell(shell_j);
    const GaussianShell &gs3 = original_bs3_->shell(shell_k);
    const GaussianShell &gs4 = original_bs4_->shell(shell_l);
    const double *xyzptr1 = gs1.center();
    double x1 = *(xyzptr1++);
    double y1 = *(xyzptr1++);
    double z1 = *(xyzptr1);
    const double *xyzptr2 = gs2.center();
    double x2 = *(xyzptr2++);
    double y2 = *(xyzptr2++);
    double z2 = *(xyzptr2);
    const double *xyzptr3 = gs3.center();
    double x3 = *(xyzptr3++);
    double y3 = *(xyzptr3++);
    double z3 = *(xyzptr3);
    const double *xyzptr4 = gs4.center();
    double x4 = *(xyzptr4++);
    double y4 = *(xyzptr4++);
    double z4 = *(xyzptr4);
    F_INT ncgto1 = 1;
    F_INT ncgto2 = 1;
    F_INT ncgto3 = 1;
    F_INT ncgto4 = 1;
    F_INT ncgto = 4;
    F_INT npgto1 = gs1.nprimitive();
    F_INT npgto2 = gs2.nprimitive();
    F_INT npgto3 = gs3.nprimitive();
    F_INT npgto4 = gs4.nprimitive();
    F_INT npgto = npgto1 + npgto2 + npgto3 + npgto4;
    F_INT am1 = gs1.am();
    F_INT am2 = gs2.am();
    F_INT am3 = gs3.am();
    F_INT am4 = gs4.am();
    ccend_[0] = npgto4;
    ccend_[1] = npgto3;
    ccend_[2] = npgto2;
    ccend_[3] = npgto1;

    // note that these are switched. shell order is reversed
    int offset_i = 0;
    int offset_j = offset_i + npgto4;
    int offset_k = offset_j + npgto3;
    int offset_l = offset_k + npgto2;

    // Copy exponents and coefficients over
    std::copy(gs4.exps(), gs4.exps()+npgto4, alpha_ + offset_i);
    std::copy(gs3.exps(), gs3.exps()+npgto3, alpha_ + offset_j);
    std::copy(gs2.exps(), gs2.exps()+npgto2, alpha_ + offset_k);
    std::copy(gs1.exps(), gs1.exps()+npgto1, alpha_ + offset_l);
    std::copy(gs4.coefs(), gs4.coefs()+npgto4, cc_ + offset_i);
    std::copy(gs3.coefs(), gs3.coefs()+npgto3, cc_ + offset_j);
    std::copy(gs2.coefs(), gs2.coefs()+npgto2, cc_ + offset_k);
    std::copy(gs1.coefs(), gs1.coefs()+npgto1, cc_ + offset_l);


#if DEBUG
    output::printf("\n\nShell (%2d %2d | %2d %2d) - center (%2d %2d | %2d %2d) - angular momentum (%d %d | %d %d)\n",
                    shell_i, shell_j, shell_k, shell_l,
                    bs1_->shell(shell_i).ncenter(), bs2_->shell(shell_j).ncenter(), bs3_->shell(shell_k).ncenter(), bs4_->shell(shell_l).ncenter(), 
                    am1, am2, am3, am4);
    output::printf("XYZ1: %16.10f %16.10f %16.10f\n", x1, y1, z1);
    output::printf("XYZ2: %16.10f %16.10f %16.10f\n", x2, y2, z2);
    output::printf("XYZ3: %16.10f %16.10f %16.10f\n", x3, y3, z3);
    output::printf("XYZ4: %16.10f %16.10f %16.10f\n", x4, y4, z4);
    output::printf("Indices  -> %d, %d\n", ccbeg_[0], ccend_[0]);

    output::printf("Number of primitives: %d %d %d %d\n", npgto1, npgto2, npgto3, npgto4);
    output::printf("Coefficients: ");
    for(int n = 0; n < npgto; ++n)
        output::printf("%14.10f ", cc_[n]);
    output::printf("\n");
    output::printf("Exponents:    ");
    for(int n = 0; n < npgto; ++n)
        output::printf("%14.10f ", alpha_[n]);
    output::printf("\n");
#endif

    F_INT nbatch = 0;
    // Call ERD.  N.B. We reverse the shell ordering, because the first index is
    // the fastest running index in the buffer, which should be l for us.
    C_ERD__GENER_ERI_BATCH(i_buffer_size_, d_buffer_size_, npgto, npgto, ncgto,
                           ncgto4, ncgto3, ncgto2, ncgto1,
                           npgto4, npgto3, npgto2, npgto1,
                           am4, am3, am2, am1,
                           x4, y4, z4, x3, y3, z3, x2, y2, z2, x1, y1, z1,
                           alpha_, cc_, ccbeg_, ccend_, spheric_, screen_,
                           iscratch_, nbatch, buffer_offset_, dscratch_);

#if DEBUG
    output::printf("Buffer offset is %d\n", buffer_offset_-1);
    output::printf("%d integrals were computed\n", nbatch);
#endif

    if(nbatch == 0){
        // The code should check the return value, and ignore the integrals in the buffer if we get here
        //::memset(target_, 0, sizeof(double)*gs1.nfunction() *gs2.nfunction() *gs3.nfunction() *gs4.nfunction());
        return 0;
    }

    if(original_bs1_->has_puream()){
        source_ = &(dscratch_[buffer_offset_-1]);
        pure_transform(shell_i, shell_j, shell_k, shell_l, 1);
    }else{
        std::copy(dscratch_ + buffer_offset_ - 1, 
                  dscratch_ + buffer_offset_ - 1 + nbatch,
                  target_);
    }
    return nbatch;

}


} // Namespace

