/*! \file
 * \brief Class for calculating libERD-based ERI (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_ERDERI_H
#define PANACHE_ERDERI_H

#include "ERDTwoElectronInt.h"

namespace panache
{


/*!
 * \brief Calculates standard two-electron ERI using libERD
 */
class ERDERI : public ERDTwoElectronInt
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
    ERDERI(const SharedBasisSet bs1,
           const SharedBasisSet bs2,
           const SharedBasisSet bs3,
           const SharedBasisSet bs4);

    /*!
     * \brief Destructor
     * 
     * Doesn't really do much
     */
    virtual ~ERDERI();
};

} // close namespace panache

#endif // PANACHE_ERDERI_H 

