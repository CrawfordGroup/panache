/*! \file
 * \brief Contnrols where PANACHE gets its two-electron integrals from
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_ERI_H
#define PANACHE_ERI_H

#include <memory>

#include "panache/Exception.h"

#ifdef PANACHE_USE_LIBINT
#include "panache/LibintERI.h"
#endif

#ifdef PANACHE_USE_LIBINT2
#include "panache/Libint2ERI.h"
#endif

#ifdef PANACHE_USE_SLOWERI
#include "panache/SlowERI.h"
#endif

#ifdef PANACHE_USE_LIBERD
#include "panache/ERDERI.h"
#endif


namespace panache {

class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;


/*!
 * \brief Obtains an ERI generator
 *
 * The generator obtained (libint1, libint2, erd, etc) depends on the options used
 * when the library was compiled.
 *
 * \param [in] bs1 Basis set on the first center
 * \param [in] bs2 Basis set on the second center
 * \param [in] bs3 Basis set on the third center
 * \param [in] bs4 Basis set on the fourth center
 * \return A generator of ERI
 */
inline SharedTwoBodyAOInt GetERI(const SharedBasisSet & bs1, const SharedBasisSet & bs2,
                                       const SharedBasisSet & bs3, const SharedBasisSet & bs4)
{
            #ifdef PANACHE_USE_LIBINT
              return SharedTwoBodyAOInt(new LibintERI(bs1, bs2, bs3, bs4));
            #endif
            #ifdef PANACHE_USE_LIBINT2
              return SharedTwoBodyAOInt(new Libint2ERI(bs1, bs2, bs3, bs4));
            #endif
            #ifdef PANACHE_USE_SLOWERI
              return SharedTwoBodyAOInt(new SlowERI(bs1, bs2, bs3, bs4));
            #endif
            #ifdef PANACHE_USE_LIBERD
              return SharedTwoBodyAOInt(new ERDERI(bs1, bs2, bs3, bs4));
            #endif
}



/*!
 * \brief Obtains an ErfERI generator
 *
 * The generator obtained (libint1, libint2, erd, etc) depends on the options used
 * when the library was compiled.
 *
 * \param [in] omega Erf \f$ \omega \f$ value
 * \param [in] bs1 Basis set on the first center
 * \param [in] bs2 Basis set on the second center
 * \param [in] bs3 Basis set on the third center
 * \param [in] bs4 Basis set on the fourth center
 * \return A generator of ErfERI
 */
inline SharedTwoBodyAOInt GetErfERI(double omega,
                                     SharedBasisSet & bs1, SharedBasisSet & bs2,
                                     SharedBasisSet & bs3, SharedBasisSet & bs4)
{
            #ifdef PANACHE_USE_LIBINT
              return SharedTwoBodyAOInt(new LibintErfERI(omega, bs1, bs2, bs3, bs4));
            #endif
            #ifdef PANACHE_USE_LIBINT2
              return SharedTwoBodyAOInt(new Libint2ErfERI(omega, bs1, bs2, bs3, bs4));
            #endif
            #ifdef PANACHE_USE_SLOWERI
              throw RuntimeError("ErfERI for SlowERI not implemented!");
            #endif
            #ifdef PANACHE_USE_LIBERD
              throw RuntimeError("ErfERI for libERD not implemented!");
            #endif
}


} // close namespace panache


#endif
