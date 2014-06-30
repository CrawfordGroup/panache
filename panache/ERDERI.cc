/*! \file
 * \brief Class for calculating libERD-based ERI (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "ERDERI.h"

namespace panache {

ERDERI::ERDERI(const SharedBasisSet bs1,
         const SharedBasisSet bs2,
         const SharedBasisSet bs3,
         const SharedBasisSet bs4)
    : ERDTwoElectronInt(bs1, bs2, bs3, bs4)
{
}

ERDERI::~ERDERI()
{
}

} // Namespace
