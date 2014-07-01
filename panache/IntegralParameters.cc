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

#include "panache/IntegralParameters.h"

namespace panache {

CorrelationFactor::CorrelationFactor(unsigned int nparam)
    : IntegralParameters(nparam)
{
}

CorrelationFactor::CorrelationFactor(const std::vector<double> & coeff, const std::vector<double> & exponent)
{
    set_params(coeff, exponent);
}

CorrelationFactor::~CorrelationFactor()
{
    delete[] coeff_;
    delete[] exponent_;
}

void CorrelationFactor::set_params(const std::vector<double> & coeff, const std::vector<double> & exponent)
{
    int nparam = coeff.size();
    if (nparam) {
        coeff_ = new double[nparam];
        exponent_ = new double[nparam];
        std::copy(coeff.begin(), coeff.end(), coeff_);
        std::copy(exponent.begin(), exponent.end(), exponent_);
    }
}

FittedSlaterCorrelationFactor::FittedSlaterCorrelationFactor(double exponent)
    : CorrelationFactor(6)
{
    // Perform the fit.
    std::vector<double> exps(6);
    std::vector<double> coeffs(6);

    slater_exponent_ = exponent;

    // The fitting coefficients
    coeffs[0] = -0.3144;
    coeffs[1] = -0.3037;
    coeffs[2] = -0.1681;
    coeffs[3] = -0.09811;
    coeffs[4] = -0.06024;
    coeffs[5] = -0.03726;

    // and the exponents
    exps[0] = 0.2209;
    exps[1] = 1.004;
    exps[2] = 3.622;
    exps[3] = 12.16;
    exps[4] = 45.87;
    exps[5] = 254.4;

    // They just need to be scaled
    double expsq = exponent * exponent;
    for(auto & it : exps)
        it *= expsq;

    set_params(coeffs, exps);
}

}
