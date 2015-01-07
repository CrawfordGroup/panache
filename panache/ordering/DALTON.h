#ifndef PANACHE_ORDERING_DALTON_H
#define PANACHE_ORDERING_DALTON_H

#include "panache/Reorder.h"

namespace panache {
namespace reorder {


class DALTON_Ordering : public Orderings
{
public:
    DALTON_Ordering(void)
    {
        // Dalton cartesian ordering seems to be the same as Psi4 ordering
        // Spherical, however, goes from -l, ..., l , where psi4 goes 0, +1, -1, ..., +l, -l
        SetSphOrder(2, {5, 3, 1, 2, 4});
        SetSphOrder(3, {7, 5, 3, 1, 2, 4, 6});
    }
};


/*!
 * \brief Returns a pointer to an Ordering object corresponding to the correct
 *        type of ordering needed
 *
 * \param [in] order Type of ordering object needed
 */
std::unique_ptr<Orderings> GetOrdering(int order);


} } // close namespace panache::reorder



#endif
