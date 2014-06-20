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

#ifndef PANACHE_LIBINTERI_H
#define PANACHE_LIBINTERI_H

#include "LibintTwoElectronInt.h"

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

class LibintERI : public LibintTwoElectronInt
{
public:
    LibintERI(const SharedBasisSet bs1,
        const SharedBasisSet bs2,
        const SharedBasisSet bs3,
        const SharedBasisSet bs4);

    virtual ~LibintERI();
};

class LibintErfERI : public LibintTwoElectronInt
{
public:
    LibintErfERI(double omega,
           const SharedBasisSet bs1,
           const SharedBasisSet bs2,
           const SharedBasisSet bs3,
           const SharedBasisSet bs4);

    virtual ~LibintErfERI();

    void setOmega(double omega);
};

}

#endif //PANACHE_LIBINTRI_H
