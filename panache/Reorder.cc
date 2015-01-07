/*! \file
 * \brief Holds information about a basis set shell (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/Reorder.h"
#include "panache/Flags.h"

// Individual ordering
#include "panache/ordering/GAMESS.h"
#include "panache/ordering/DALTON.h"


namespace panache {
namespace reorder {

std::unique_ptr<Orderings> GetOrdering(int order)
{
    switch(order)
    {
        case BSORDER_PSI4:
            return std::unique_ptr<Orderings>(); // null
            break;
        case BSORDER_GAMESS:
            return std::unique_ptr<Orderings>(new reorder::GAMESS_Ordering());
            break;
        case BSORDER_DALTON:
            return std::unique_ptr<Orderings>(new reorder::DALTON_Ordering());
            break;
        default:
            throw RuntimeError("Error - unknown ordering given to GetOrdering!");
    }
}


std::unique_ptr<CNorm> GetCNorm(int order)
{
    switch(order)
    {
        case BSORDER_PSI4:
            return std::unique_ptr<CNorm>(); // null
            break;
        case BSORDER_GAMESS:
            return std::unique_ptr<CNorm>(new reorder::GAMESS_CNorm());
            break;
        case BSORDER_DALTON:
            return std::unique_ptr<CNorm>(); // null
            break;
        default:
            throw RuntimeError("Error - unknown ordering given to GetCNorm!");
    }
}



}} // close namespace panache::reorder
