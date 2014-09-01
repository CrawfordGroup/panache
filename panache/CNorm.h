#ifndef PANACHE_CNORM_H
#define PANACHE_CNORM_H

#include <vector>
#include <array>
#include <functional>
#include <sstream>

#ifndef MAX_CNORM_AM
#define MAX_CNORM_AM 10
#endif

using std::vector;
using std::array;
using std::function;


namespace panache {
namespace reorder {


class CNorm
{
private:
    // we keep a vector so we can have a homogeneous container
    array<vector<double>, MAX_CNORM_AM > cart_cnorm_;
    array<vector<double>, MAX_CNORM_AM > sph_cnorm_;

    void CheckAM(int am, const char * func) const
    {
        std::stringstream ss;
        ss << "Cannot access C normalization for am = " << am << " in function " << func << " - MAX is " << MAX_CNORM_AM << "\n";
        if(am > MAX_CNORM_AM)
            throw RuntimeError(ss.str());
    }

public:
    void SetCartCNorm(int am, vector<double> order)
    {
        CheckAM(am, __FUNCTION__);

        cart_cnorm_[am].assign(order.begin(), order.end());
    }

    void SetSphCNorm(int am, vector<double> order)
    {
        CheckAM(am, __FUNCTION__);

        sph_cnorm_[am].assign(order.begin(), order.end());
    }


    vector<double> GetCartCNorm(int am) const
    {
        CheckAM(am, __FUNCTION__);

        return cart_cnorm_[am];
    }

    vector<double> GetSphCNorm(int am) const
    {
        CheckAM(am, __FUNCTION__);

        return sph_cnorm_[am];
    }

    vector<double> GetCNorm(bool ispure, int am) const
    {
        CheckAM(am, __FUNCTION__);

        return ispure ? sph_cnorm_[am] : cart_cnorm_[am];
    }


    bool NeedsCartCNorm(int am) const
    {
        CheckAM(am, __FUNCTION__);

        return cart_cnorm_[am].size();
    }

    bool NeedsSphCNorm(int am) const
    {
        CheckAM(am, __FUNCTION__);

        return sph_cnorm_[am].size();
    }

    bool NeedsCNorm(bool ispure, int am) const
    {
        CheckAM(am, __FUNCTION__);

        return ispure ? sph_cnorm_[am].size() : cart_cnorm_[am].size();
    }
};



class GAMESS_CNorm : public CNorm
{
public:
    GAMESS_CNorm(void)
    {
        const double sqrt3 = sqrt(3.0);

        // don't ask
        SetCartCNorm(2, {1.0, 1.0, 1.0, sqrt3, sqrt3, sqrt3} );
    }
};



} } // close namespace panache::reorder



#endif
