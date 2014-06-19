/*! \file
 * \brief Defines some integer types to 32- or 64-bit
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_INT_T_H
#define PANACHE_INT_T_H

#include <cstdint>

/*! \def int_t
 *
 *  \brief Is set to either int64_t or int32_t, depending on the USE_64PANACHE option.
 *         See \ref compiling_f64bit_sec
 */

#ifdef USE_64PANACHE
    #define int_t int64_t
#else
    #define int_t int32_t
#endif


/*! \def lapack_int_t
 *
 *  \brief Is set to either int64_t or int32_t, depending on the USE_64LAPACK option.
 *         See \ref compiling_lapack64
 */
#ifdef USE_64LAPACK
  #define lapack_int_t int64_t
#else
  #define lapack_int_t int32_t
#endif



#endif
