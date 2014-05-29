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

#ifndef PANACHE_LIBINT2ERI_H
#define PANACHE_LIBINT2ERI_H

#include "Libint2TwoElectronInt.h"

namespace panache
{

class BasisSet;
class GaussianShell;
class TwoBodyAOInt;
class IntegralFactory;
class SphericalTransform;
class Fjt;
class AOShellCombinationsIterator;
class CorrelationFactor;

class Libint2ERI : public Libint2TwoElectronInt
{
public:
    Libint2ERI(const std::shared_ptr<BasisSet> bs1,
        const std::shared_ptr<BasisSet> bs2,
        const std::shared_ptr<BasisSet> bs3,
        const std::shared_ptr<BasisSet> bs4);

    virtual ~Libint2ERI();
};

class Libint2ErfERI : public Libint2TwoElectronInt
{
public:
    Libint2ErfERI(double omega,
           const std::shared_ptr<BasisSet> bs1,
           const std::shared_ptr<BasisSet> bs2,
           const std::shared_ptr<BasisSet> bs3,
           const std::shared_ptr<BasisSet> bs4);

    virtual ~Libint2ErfERI();

    void setOmega(double omega);
};

}

#endif //PANACHE_LIBINT2ERI_H
