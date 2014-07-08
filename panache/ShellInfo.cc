/*! \file
 * \brief Holds information about a basis set shell (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <cstdlib>
#include <cmath>

#include "panache/ShellInfo.h"
#include "panache/BasisFunctionMacros.h" // for INT_NCART and INT_NFUNC
#include "panache/Math.h"
#include "panache/PhysConst.h"
#include "panache/Output.h"

namespace panache {

ShellInfo::ShellInfo(int am, const std::vector<double> &c,
                             const std::vector<double> &e, GaussianType pure,
                             int nc, const Vector3 &center,
                             PrimitiveType pt)
    : l_(am), puream_(pure), exp_(e), coef_(c), original_coef_(c),
      nc_(nc), center_(center)
{
    ncartesian_ = INT_NCART(l_);
    nfunction_  = INT_NFUNC(puream_, l_);

    // Compute the normalization constants
    if (pt == Unnormalized){
        normalize_shell();
    }

    // by default, coef_ = original_coef_ = c

}


double ShellInfo::primitive_normalization(int p)
{
    double tmp1 = l_ + 1.5;
    double g = 2.0 * exp_[p];
    double z = pow(g, tmp1);
    double normg = sqrt( (pow(2.0, l_) * z) / (M_PI * sqrt(M_PI) * math::double_factorial_nminus1(2*l_)));
    return normg;
}


void ShellInfo::normalize_shell()
{
    double m = (double)l_+1.5;
    double sum = 0.0;
    for(int j = 0; j < nprimitive(); j++){
        for(int k = 0; k <= j; k++){
            double a1 = exp_[j];
            double a2 = exp_[k];
            double temp = (original_coef(j) * original_coef(k));
            double temp2 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
            temp2 = pow(temp2, m);
            temp = temp * temp2;
            sum = sum + temp;
            if(j != k)
                sum = sum + temp;
        }
    }

#if PANACHE_USE_LIBERD
    double prefac = 1.0;
    if(l_ > 1)
        prefac = pow(2.0, 2*l_) / math::double_factorial_nminus1(2*l_);
    double norm = sqrt(prefac / sum);
    for(int i = 0; i < nprimitive(); i++){
        coef_[i] *= norm * pow(exp_[i], 0.5*m);
    }
#else
    for (int i = 0; i < nprimitive(); ++i) {
        double norm = sqrt(1.0 / sum);
        coef_[i] *= norm * primitive_normalization(i);
    }
#endif
}

int ShellInfo::nfunction() const
{
    return INT_NFUNC(puream_, l_);
}

int ShellInfo::nprimitive() const
{
    return exp_.size();
}

const Vector3& ShellInfo::center() const
{
    return center_;
}

void ShellInfo::print(void) const
{
    output::printf("    %c %3d 1.00\n", AMCHAR(), nprimitive());
    for (int K = 0; K < nprimitive(); K++)
        output::printf("               %20.8f %20.8f\n",exp_[K], original_coef_[K]);
}

const char *ShellInfo::amtypes = "spdfghiklmnopqrtuvwxyz";
const char *ShellInfo::AMTYPES = "SPDFGHIKLMNOPQRTUVWXYZ";

} // close namespace panache
