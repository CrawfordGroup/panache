#ifndef PANACHE_REORDER_H
#define PANACHE_REORDER_H

#include "panache/Orderings.h"
#include "panache/CNorm.h"


namespace panache {
namespace reorder {

/*!
 * \brief Returns a pointer to an Ordering object corresponding to the correct
 *        type of ordering needed
 *
 * \param [in] order Type of ordering object needed
 */
std::unique_ptr<Orderings> GetOrdering(int order);


/*!
 * \brief Returns a pointer to an CNorm object corresponding to the correct
 *        type of c-matrix renormalization needed
 *
 * \param [in] order Type of cnorm object needed
 */
std::unique_ptr<CNorm> GetCNorm(int order);


} } // close namespace panache::reorder



#endif
