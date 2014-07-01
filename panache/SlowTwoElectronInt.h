/*! \file
 * \brief Base class for SlowERI-based ERI (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_SLOWTWOELECTRONINT_H
#define PANACHE_SLOWTWOELECTRONINT_H

#include "panache/TwoBodyAOInt.h"
#include "panache/SlowERIBase.h"

namespace panache
{


/*!
 * \brief Implements the SlowERI backend for calculation
 *        of electron repulsion integrals.
 */
class SlowTwoElectronInt : public TwoBodyAOInt
{
private:
    SlowERIBase sloweri_; //!< SlowERIBase object that actually does the calculating


    /*!
     * \brief Computes the ERIs between four shells.
     *
     * \param [in] sh1 Shell 1 of the quartet
     * \param [in] sh2 Shell 2 of the quartet
     * \param [in] sh3 Shell 3 of the quartet
     * \param [in] sh4 Shell 4 of the quartet
     * \return Number of computed integrals
     */
    size_t compute_quartet(int sh1, int sh2, int sh3, int sh4);


public:
    /*!
     * \brief Constructor 
     *
     * \param [in] bs1 The basis set on the 1st center
     * \param [in] bs2 The basis set on the 2nd center
     * \param [in] bs3 The basis set on the 3rd center
     * \param [in] bs4 The basis set on the 4th center
     */
    SlowTwoElectronInt(const SharedBasisSet bs1,
                       const SharedBasisSet bs2,
                       const SharedBasisSet bs3,
                       const SharedBasisSet bs4);

    virtual ~SlowTwoElectronInt();

    // See TwoBodyAOInt::compute_shell
    virtual size_t compute_shell(const AOShellCombinationsIterator& shellIter);

    // See TwoBodyAOInt::compute_shell
    virtual size_t compute_shell(int sh1, int sh2, int sh3, int sh4);
};

}

#endif //PANACHE_SLOWTWOELECTRONINT_H


