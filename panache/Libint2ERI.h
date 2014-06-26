/*! \file
 * \brief Class for calculating libint2-based ERI (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_LIBINT2ERI_H
#define PANACHE_LIBINT2ERI_H

#include "Libint2TwoElectronInt.h"

namespace panache
{


/*!
 * \brief Calculates standard two-electron ERI using libint2
 */
class Libint2ERI : public Libint2TwoElectronInt
{
public:
    /*!
     * \brief Constructor 
     *
     * Also initializes the libint2 library and allocates memory
     *
     * \param [in] bs1 The basis set on the 1st center
     * \param [in] bs2 The basis set on the 2nd center
     * \param [in] bs3 The basis set on the 3rd center
     * \param [in] bs4 The basis set on the 4th center
     */
    Libint2ERI(const SharedBasisSet bs1,
               const SharedBasisSet bs2,
               const SharedBasisSet bs3,
               const SharedBasisSet bs4);

    /*!
     * \brief Destructor
     * 
     * Frees memory and closes the libint2 library
     */
    virtual ~Libint2ERI();
};



/*!
 * \brief Calculates attenuated two-electron ERI using libint
 */
class Libint2ErfERI : public Libint2TwoElectronInt
{
public:

    /*!
     * \brief Constructor 
     *
     * Also initializes the libint2 library and allocates memory
     *
     * \param [in] omega Attenuation parameter
     * \param [in] bs1 The basis set on the 1st center
     * \param [in] bs2 The basis set on the 2nd center
     * \param [in] bs3 The basis set on the 3rd center
     * \param [in] bs4 The basis set on the 4th center
     */
    Libint2ErfERI(double omega,
           const SharedBasisSet bs1,
           const SharedBasisSet bs2,
           const SharedBasisSet bs3,
           const SharedBasisSet bs4);


    /*!
     * \brief Destructor
     * 
     * Frees memory and destroys the libint object
     */
    virtual ~Libint2ErfERI();


    /*!
     * \brief Changes the attenuation parameter omega
     *
     * \param [in] omega The new omega to change to
     */ 
    void setOmega(double omega);
};

}

#endif //PANACHE_LIBINT2ERI_H
