/*! \file
 * \brief Class for calculating SlowERI-based ERI (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_SLOWERI_H
#define PANACHE_SLOWERI_H

#include "SlowTwoElectronInt.h"

namespace panache
{


/*!
 * \brief Calculates standard two-electron ERI using libint
 */
class SlowERI : public SlowTwoElectronInt
{
public:
    /*!
     * \brief Constructor 
     *
     * \param [in] bs1 The basis set on the 1st center
     * \param [in] bs2 The basis set on the 2nd center
     * \param [in] bs3 The basis set on the 3rd center
     * \param [in] bs4 The basis set on the 4th center
     */
    SlowERI(const SharedBasisSet bs1,
           const SharedBasisSet bs2,
           const SharedBasisSet bs3,
           const SharedBasisSet bs4);

    /*!
     * \brief Destructor
     *
     * Doesn't really do anything
     */ 
    virtual ~SlowERI();
};

} // close namespace panache

#endif // PANACHE_SLOWERI_H 

