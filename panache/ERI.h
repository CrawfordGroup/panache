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



using std::shared_ptr;

namespace panache {


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
inline shared_ptr<TwoBodyAOInt> GetERI(shared_ptr<BasisSet> & bs1, shared_ptr<BasisSet> & bs2,
                                       shared_ptr<BasisSet> & bs3, shared_ptr<BasisSet> & bs4)
{
            #ifdef PANACHE_USE_LIBINT
              return std::shared_ptr<TwoBodyAOInt>(new LibintERI(bs1, bs2, bs3, bs4));
            #endif
            #ifdef PANACHE_USE_LIBINT2
              return std::shared_ptr<TwoBodyAOInt>(new Libint2ERI(bs1, bs2, bs3, bs4));
            #endif
            #ifdef PANACHE_USE_SLOWERI
              return std::shared_ptr<TwoBodyAOInt>(new SlowERI(bs1, bs2, bs3, bs4));
            #endif
            #ifdef PANACHE_USE_LIBERD
              return std::shared_ptr<TwoBodyAOInt>(new ERDERI(bs1, bs2, bs3, bs4));
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
inline shared_ptr<TwoBodyAOInt> GetErfERI(double omega,
                                          shared_ptr<BasisSet> & bs1, shared_ptr<BasisSet> & bs2,
                                          shared_ptr<BasisSet> & bs3, shared_ptr<BasisSet> & bs4)
{
            #ifdef PANACHE_USE_LIBINT
              return std::shared_ptr<TwoBodyAOInt>(new LibintErfERI(omega, bs1, bs2, bs3, bs4));
            #endif
            #ifdef PANACHE_USE_LIBINT2
              return std::shared_ptr<TwoBodyAOInt>(new Libint2ErfERI(omega, bs1, bs2, bs3, bs4));
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
