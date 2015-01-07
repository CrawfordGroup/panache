#ifndef PANACHE_ORDERINGS_H
#define PANACHE_ORDERINGS_H

#include <vector>
#include <array>
#include <sstream>
#include <memory>

#include "panache/BasisFunctionMacros.h"
#include "panache/Exception.h"

#ifndef MAX_REORDER_AM
#define MAX_REORDER_AM 10
#endif

using std::vector;
using std::array;


namespace panache {
namespace reorder {


struct PointerMap
{
    int start;
    vector<unsigned short> order;

    PointerMap(int start, const vector<unsigned short> & order)
        : start(start), order(order)
    { }

};



class Orderings
{
private:
    // we keep a vector so we can have a homogeneous container
    array<vector<unsigned short>, MAX_REORDER_AM+1 > cart_orderings_;
    array<vector<unsigned short>, MAX_REORDER_AM+1 > sph_orderings_;
    array<vector<unsigned short>, MAX_REORDER_AM+1 > inv_cart_orderings_;
    array<vector<unsigned short>, MAX_REORDER_AM+1 > inv_sph_orderings_;

    void CheckAM(int am, const char * func) const
    {
        std::stringstream ss;
        ss << "Cannot access ordering for am = " << am << " in function " << func << " - MAX is " << MAX_REORDER_AM << "\n";
        if(am > MAX_REORDER_AM)
            throw RuntimeError(ss.str());
    }

public:
    void SetCartOrder(int am, const vector<unsigned short> & order)
    {
        CheckAM(am, __FUNCTION__);

        cart_orderings_[am].assign(order.begin(), order.end());

        inv_cart_orderings_[am].resize(order.size());
        for(size_t i = 0; i < order.size(); i++)
            inv_cart_orderings_[am][order[i]-1] = i+1;
    }

    void SetSphOrder(int am, const vector<unsigned short> & order)
    {
        CheckAM(am, __FUNCTION__);

        sph_orderings_[am].assign(order.begin(), order.end());

        inv_sph_orderings_[am].resize(order.size());
        for(size_t i = 0; i < order.size(); i++)
            inv_sph_orderings_[am][order[i]-1] = i+1;
    }


    vector<unsigned short> GetCartOrder(int am) const
    {
        CheckAM(am, __FUNCTION__);
        return cart_orderings_[am];
    }

    vector<unsigned short> GetSphOrder(int am) const
    {
        CheckAM(am, __FUNCTION__);
        return sph_orderings_[am];
    }

    vector<unsigned short> GetOrder(bool ispure, int am) const
    {
        CheckAM(am, __FUNCTION__);
        return ispure ? sph_orderings_[am] : cart_orderings_[am];
    }

    vector<unsigned short> GetInvCartOrder(int am) const
    {
        CheckAM(am, __FUNCTION__);
        return inv_cart_orderings_[am];
    }

    vector<unsigned short> GetInvSphOrder(int am) const
    {
        CheckAM(am, __FUNCTION__);
        return inv_sph_orderings_[am];
    }

    vector<unsigned short> GetInvOrder(bool ispure, int am) const
    {
        CheckAM(am, __FUNCTION__);
        return ispure ? inv_sph_orderings_[am] : inv_cart_orderings_[am];
    }


    bool NeedsCartReordering(int am) const
    {
        CheckAM(am, __FUNCTION__);
        return cart_orderings_[am].size();
    }

    bool NeedsSphReordering(int am) const
    {
        CheckAM(am, __FUNCTION__);
        return sph_orderings_[am].size();
    }

    bool NeedsReordering(bool ispure, int am) const
    {
        CheckAM(am, __FUNCTION__);
        return ispure ? sph_orderings_[am].size() : cart_orderings_[am].size();
    }

    bool NeedsInvCartReordering(int am) const
    {
        CheckAM(am, __FUNCTION__);
        return inv_cart_orderings_[am].size();
    }

    bool NeedsInvSphReordering(int am) const
    {
        CheckAM(am, __FUNCTION__);
        return inv_sph_orderings_[am].size();
    }

    bool NeedsInvReordering(bool ispure, int am) const
    {
        CheckAM(am, __FUNCTION__);
        return ispure ? inv_sph_orderings_[am].size() : inv_cart_orderings_[am].size();
    }

};


} } // close namespace panache::reorder



#endif
