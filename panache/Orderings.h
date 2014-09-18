#ifndef PANACHE_REORDER_H
#define PANACHE_REORDER_H

#include <vector>
#include <array>
#include <functional>
#include <sstream>

#include "BasisFunctionMacros.h"

#ifndef MAX_REORDER_AM
#define MAX_REORDER_AM 10
#endif

using std::vector;
using std::array;
using std::function;


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



class GAMESS_Ordering : public Orderings
{
public:
    GAMESS_Ordering(void)
    {
        // don't ask
        SetCartOrder(2, {1, 4, 6, 2, 3, 5} );
        SetCartOrder(3, {1, 7, 10, 2, 3, 4, 8, 6, 9, 5} );
        SetCartOrder(4, {1, 11, 15, 2, 3, 7, 12, 10, 14, 4, 6, 13, 5, 8, 9} );
        SetCartOrder(5, {1, 16, 21, 2, 3, 11, 17, 15, 20, 4, 6, 7, 18, 10, 19, 8, 9, 13, 9, 13, 9} );
        SetCartOrder(6, {1, 22, 28, 2, 3, 16, 23, 21, 27, 4, 6, 11, 24, 15, 26, 7, 10, 25, 8, 9, 12, 18, 14, 19, 13, 13, 13, 13} );
        SetCartOrder(7, {1, 29, 36, 2, 3, 22, 30, 28, 35, 4, 6, 16, 31, 21, 34, 7, 10, 11, 32, 15, 33, 12, 14, 25, 13, 18, 19, 18, 19, 18, 19, 18, 19, 18, 19, 18} );
        SetCartOrder(8, {1, 37, 45, 2, 3, 29, 38, 36, 44, 4, 6, 22, 39, 28, 43, 7, 10, 16, 40, 21, 42, 11, 15, 41, 12, 14, 17, 32, 20, 33, 18, 19, 25, 19, 25, 19, 25, 19, 25, 19, 25, 19, 25, 19, 25} );

        for(int am = 1; am <= MAX_REORDER_AM; am++)
        {
            int nfunc = INT_NPURE(am);
            vector<unsigned short> re(nfunc);

            re[0] = am+1;
            int j = 0;
            for(int i = 1; i <= (nfunc-1)/2; i++)
            {
                re[++j] = re[0] + i;
                re[++j] = re[0] - i;
            }

            SetSphOrder(am, re);
        }

    }
};



} } // close namespace panache::reorder



#endif
