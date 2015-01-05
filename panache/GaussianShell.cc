/*! \file
 * \brief Holds information about a basis set shell (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <cstdlib>

#include "panache/GaussianShell.h"
#include "panache/BasisFunctionMacros.h"
#include "panache/Output.h"

namespace panache
{

GaussianShell::GaussianShell(int am, int nprimitive, const double *oc,
                             const double *c,
                             const double *e, ShellInfo::GaussianType pure,
                             int nc, const double *center, int start)
    : l_(am), puream_(pure), 
      exp_(e),
      original_coef_(oc),
      coef_(c), 
      nc_(nc), center_(center), start_(start), nprimitive_(nprimitive)
{
    ncartesian_ = INT_NCART(l_);
    nfunction_  = INT_NFUNC(puream_, l_);
}

int GaussianShell::nfunction() const
{
    return INT_NFUNC(puream_, l_);
}

int GaussianShell::nprimitive() const
{
    return nprimitive_;
}

const double* GaussianShell::center() const
{
    return center_;
}

void GaussianShell::print(void) const
{
    output::printf("    %c %3d   %s\n", AMCHAR(), nprimitive(), is_pure() ? "sph" : "cart");
    for (int K = 0; K < nprimitive(); K++)
        output::printf("               %20.8f %20.8f\n",exp_[K], coef_[K]);
}

const char *GaussianShell::amtypes = "spdfghiklmnopqrtuvwxyz";
const char *GaussianShell::AMTYPES = "SPDFGHIKLMNOPQRTUVWXYZ";

} // close namespace panache


