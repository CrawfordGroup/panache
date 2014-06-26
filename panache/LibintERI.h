/*! \file
 * \brief Class for calculating libint1-based ERI (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_LIBINTERI_H
#define PANACHE_LIBINTERI_H

#include "LibintTwoElectronInt.h"

namespace panache
{


/*!
 * \brief Calculates standard two-electron ERI using libint
 */
class LibintERI : public LibintTwoElectronInt
{
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
    LibintERI(const SharedBasisSet bs1,
              const SharedBasisSet bs2,
              const SharedBasisSet bs3,
              const SharedBasisSet bs4);

    /*!
     * \brief Destructor
     * 
     * Frees memory and destroys the libint object
     */
    virtual ~LibintERI();
};



/*!
 * \brief Calculates attenuated two-electron ERI using libint
 */
class LibintErfERI : public LibintTwoElectronInt
{
public:

    /*!
     * \brief Constructor 
     *
     * Also initializes the libint library and allocates memory
     *
     * \param [in] omega Attenuation parameter
     * \param [in] bs1 The basis set on the 1st center
     * \param [in] bs2 The basis set on the 2nd center
     * \param [in] bs3 The basis set on the 3rd center
     * \param [in] bs4 The basis set on the 4th center
     */
    LibintErfERI(double omega,
                 const SharedBasisSet bs1,
                 const SharedBasisSet bs2,
                 const SharedBasisSet bs3,
                 const SharedBasisSet bs4);


    /*!
     * \brief Destructor
     * 
     * Frees memory and destroys the libint object
     */
    virtual ~LibintErfERI();


    /*!
     * \brief Changes the attenuation parameter omega
     *
     * \param [in] omega The new omega to change to
     */ 
    void setOmega(double omega);
};

}

#endif //PANACHE_LIBINTRI_H
