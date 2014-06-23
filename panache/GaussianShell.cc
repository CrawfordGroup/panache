/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <cstdlib>

#include "GaussianShell.h"
#include "BasisFunctionMacros.h"
#include "Output.h"

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
    output::printf("    %c %3d 1.00\n", AMCHAR(), nprimitive());
    for (int K = 0; K < nprimitive(); K++)
        output::printf("               %20.8f %20.8f\n",exp_[K], coef_[K]);
}

const char *GaussianShell::amtypes = "spdfghiklmnopqrtuvwxyz";
const char *GaussianShell::AMTYPES = "SPDFGHIKLMNOPQRTUVWXYZ";

} // close namespace panache


