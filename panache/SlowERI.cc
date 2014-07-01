/*! \file
 * \brief Class for calculating SlowERI-based ERI (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/SlowERI.h"

namespace panache {

SlowERI::SlowERI(const SharedBasisSet bs1,
         const SharedBasisSet bs2,
         const SharedBasisSet bs3,
         const SharedBasisSet bs4)
    : SlowTwoElectronInt(bs1, bs2, bs3, bs4)
{
// nothing really needs to be done
}

SlowERI::~SlowERI()
{
}

} // Namespace
