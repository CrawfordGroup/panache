/*! \file
 * \brief Defines int_t to 32- or 64-bit
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_INT_T_H
#define PANACHE_INT_T_H

#include <cstdint>

/*! \def int_t
 *
 *  \brief Set to either int64_t or int32_t, depending on the USE_64PANACHE option
 *
 */

#ifdef USE_64PANACHE
    #define int_t int64_t
#else
    #define int_t int32_t
#endif

#endif
