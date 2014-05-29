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

#include "LibintERI.h"
#include "Fjt.h"
#include "BasisSet.h"

namespace panache
{

/////////
// Normal two-electron repulsion integrals
/////////

ERI::ERI(const std::shared_ptr<BasisSet> bs1,
         const std::shared_ptr<BasisSet> bs2,
         const std::shared_ptr<BasisSet> bs3,
         const std::shared_ptr<BasisSet> bs4)
    : TwoElectronInt(bs1, bs2, bs3, bs4)
{
    fjt_ = new Taylor_Fjt(basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          1, 1e-15);
}

ERI::~ERI()
{
    delete fjt_;
}

/////////
// ErfERI
/////////

ErfERI::ErfERI(double omega,
               const std::shared_ptr<BasisSet> bs1,
               const std::shared_ptr<BasisSet> bs2,
               const std::shared_ptr<BasisSet> bs3,
               const std::shared_ptr<BasisSet> bs4)
    : TwoElectronInt(bs1, bs2, bs3, bs4)
{
    fjt_ = new ErfFundamental(omega,
                              basis1()->max_am() +
                              basis2()->max_am() +
                              basis3()->max_am() +
                              basis4()->max_am() +
                              1);
}

ErfERI::~ErfERI()
{
    delete fjt_;
}

void ErfERI::setOmega(double omega)
{
    (static_cast<ErfFundamental*>(fjt_))->setOmega(omega);
}

} // close namespace panache
