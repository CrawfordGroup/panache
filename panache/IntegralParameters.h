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

#ifndef PANACHE_INTEGRALPARAMETERS_H
#define PANACHE_INTEGRALPARAMETERS_H

#include <memory>
#include <vector>

namespace panache {

class IntegralParameters
{
private:
    unsigned int nparam_;
public:
    IntegralParameters(unsigned int nparam=0)
        : nparam_(nparam) { }
    virtual ~IntegralParameters() {}

    unsigned int nparam() const { return nparam_; }
};

class CorrelationFactor : public IntegralParameters
{
private:
    double *coeff_;
    double *exponent_;

public:
    CorrelationFactor(unsigned int nparam);
    CorrelationFactor(const std::vector<double> & coeff,
                      const std::vector<double> & exponent);
    virtual ~CorrelationFactor();

    void set_params(const std::vector<double> & coeff,
                    const std::vector<double> & exponent);
    double *exponent() const { return exponent_; }
    double *coeff() const { return coeff_; }
};

typedef std::shared_ptr<CorrelationFactor> SharedCorrelationFactor;

class FittedSlaterCorrelationFactor : public CorrelationFactor
{
private:
    double slater_exponent_;

public:
    FittedSlaterCorrelationFactor(double exponent);
    double exponent(){return slater_exponent_;}
};

}

#endif //PANACHE_INTEGRALPARAMETERS_H
