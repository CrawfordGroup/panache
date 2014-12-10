/*! \file
 * \brief Some basis math functions (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_MATH_H
#define PANACHE_MATH_H

#define FACTORIAL_TABLE_SIZE 15
#define DBL_FACTORIAL_TABLE_SIZE 15
#define PACKED_ROW_TABLE_SIZE 2000

#include <utility>
#include <cstdint>

namespace panache {

/*!
 * \brief Functions for some simple math
 */
namespace math {

/*!
 * \brief Initialize some lookup tables
 */
void initmath(void);


/*!
 * \brief Simple factorial
 *
 * Integer \p n must be positive!
 *
 * \param [in] n Integer to get the factorial of
 * \return A floating point (double) of n!
 */
double factorial(int n);


/*!
 * \brief Simple double factorial
 *
 * Integer \p n must be positive!
 *
 * \param [in] n Integer to get the double factorial of
 * \return A floating point (double) of n!!
 */
double double_factorial(int n);



/*!
 * \brief Simple double factorial of (\p n - 1)
 *
 * Integer \p n must be positive!
 *
 * \param [in] n Integer to get the double factorial of
 * \return (n-1)!!
 */
double double_factorial_nminus1(int n);



/*!
 * \brief Combinations function
 *
 *  Calculates \f$ _nC_k \f$
 *
 *
 * \param [in] n Total number
 * \param [in] k Choose k
 * \return \f$ _nC_k \f$
 */
double combinations(int n, int k);



/*!
 * \brief Decompose a packed index ij into i and j
 *
 * \param [in] ij The packed index
 * \return separate i and j values
 */
std::pair<int64_t, int64_t> decomposeij_packed(int64_t ij); 

}} // close namespace panache::Math

#endif //PANACHE_MATH_H
