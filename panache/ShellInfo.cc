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


void ShellInfo::contraction_normalization()
{
    int i, j;
    double e_sum = 0.0, g, z;

    for (i=0; i<nprimitive(); ++i) {
        for (j=0; j<nprimitive(); ++j) {
            g = exp_[i] + exp_[j];
            z = pow(g, l_+1.5);
            e_sum += coef_[i] * coef_[j] / z;
        }
    }

    double tmp = ((2.0*M_PI/M_2_SQRTPI) * math::double_factorial_nminus1(2*l_))/pow(2.0, l_);
    double norm = sqrt(1.0 / (tmp*e_sum));

    // Set the normalization
    for (i=0; i<nprimitive(); ++i)
        coef_[i] *= norm;
}

void ShellInfo::normalize_shell()
{
#if PANACHE_USE_LIBERD
    double sum = 0.0;
    for(int j = 0; j < nprimitive(); j++){
        for(int k = 0; k <= j; k++){
            double a1 = exp_[j];
            double a2 = exp_[k];
            double temp = (original_coef(j) * original_coef(k));
            double temp2 = ((double) l_ + 1.5);
            double temp3 = (2.0 * sqrt(a1 * a2) / (a1 + a2));
            temp3 = pow(temp3, temp2);
            temp = temp * temp3;
            sum = sum + temp;
            if(j != k)
                sum = sum + temp;
        }
    }
    double prefac = 1.0;
    if(l_ > 1)
        prefac = pow(2.0, 2*l_) / math::double_factorial_nminus1(2*l_);
    double norm = sqrt(prefac / sum);
    for(int j = 0; j < nprimitive(); j++){
        coef_[j] *= norm;
    }
#else
    for (int i = 0; i < nprimitive(); ++i) {
        double normalization = primitive_normalization(i);
        coef_[i] *= normalization;
    }
    contraction_normalization();
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
