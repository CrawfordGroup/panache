/*! \file
 * \brief Base class for libint1-based ERI (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_LIBINTTWOELECTRONINT_H
#define PANACHE_LIBINTTWOELECTRONINT_H

#include <libint/libint.h>
#include "TwoBodyAOInt.h"

namespace panache {


class Fjt;
class AOShellCombinationsIterator;


/*!
 * \brief Implements the Libint1 backend for calculation
 *        of electron repulsion integrals.
 */
class LibintTwoElectronInt : public TwoBodyAOInt
{
protected:
    
    Libint_t libint_;  //!< Libint object.

    int max_cart_;  //!< Maximum cartesian class size.
    
    Fjt *fjt_;  //!< Computes the fundamental


    int osh1_,  //!< Original shell 1 index requested
        osh2_,  //!< Original shell 2 index requested
        osh3_,  //!< Original shell 3 index requested
        osh4_;  //!< Original shell 4 index requested

    bool p13p24_, //!< Indices 1 & 3 and 2 & 4 were permuted
         p12_,    //!< Indices 1 and 2 were permuted
         p34_;    //!< Indices 3 and 4 were permuted



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
     * Also initializes the libint library and allocates memory
     *
     * \param [in] bs1 The basis set on the 1st center
     * \param [in] bs2 The basis set on the 2nd center
     * \param [in] bs3 The basis set on the 3rd center
     * \param [in] bs4 The basis set on the 4th center
     */
    LibintTwoElectronInt(const SharedBasisSet bs1,
                         const SharedBasisSet bs2,
                         const SharedBasisSet bs3,
                         const SharedBasisSet bs4);


    virtual ~LibintTwoElectronInt();


    // See TwoBodyAOInt::compute_shell
    virtual size_t compute_shell(const AOShellCombinationsIterator& shellIter);

    // See TwoBodyAOInt::compute_shell
    virtual size_t compute_shell(int sh1, int sh2, int sh3, int sh4);

};

}

#endif //PANACHE_LIBINTTWOELECTRONINT_H
