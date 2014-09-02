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
        // don't ask
        SetCartCNorm(2,
                        {1,
                        1,
                        1,
                        1.73205080756888,
                        1.73205080756888,
                        1.73205080756888
                        }
                    ); 

        SetCartCNorm(3,
                        {1,
						 1,
						 1,
						 2.23606797749979,
						 2.23606797749979,
						 2.23606797749979,
						 2.23606797749979,
						 2.23606797749979,
						 2.23606797749979,
						 3.87298334620742
						}
                    );
        SetCartCNorm(4,
						{1,
						 1,
						 1,
						 2.64575131106459,
						 2.64575131106459,
						 2.64575131106459,
						 2.64575131106459,
						 2.64575131106459,
						 2.64575131106459,
						 3.41565025531987,
						 3.41565025531987,
						 3.41565025531987,
						 5.91607978309962,
						 5.91607978309962,
						 5.91607978309962
						}
                    );

        SetCartCNorm(5,
						{1,
						 1,
						 1,
						 3,
						 3,
						 3,
						 3,
						 3,
						 3,
						 4.58257569495584,
						 4.58257569495584,
						 4.58257569495584,
						 4.58257569495584,
						 4.58257569495584,
						 4.58257569495584,
						 10.2469507659596,
						 10.2469507659596,
						 10.2469507659596,
						 10.2469507659596,
						 10.2469507659596,
						 10.2469507659596
						}
                    );

        SetCartCNorm(6,
						{1,
						 1,
						 1,
						 3.3166247903554,
						 3.3166247903554,
						 3.3166247903554,
						 3.3166247903554,
						 3.3166247903554,
						 3.3166247903554,
						 5.74456264653803,
						 5.74456264653803,
						 5.74456264653803,
						 5.74456264653803,
						 5.74456264653803,
						 5.74456264653803,
						 6.79705818718657,
						 6.79705818718657,
						 6.79705818718657,
						 15.1986841535707,
						 15.1986841535707,
						 15.1986841535707,
						 15.1986841535707,
						 15.1986841535707,
						 15.1986841535707,
						 19.6214168703486,
						 19.6214168703486,
						 19.6214168703486,
						 19.6214168703486
						}
                    );

        SetCartCNorm(7,
						{1,
						 1,
						 1,
						 3.60555127546399,
						 3.60555127546399,
						 3.60555127546399,
						 3.60555127546399,
						 3.60555127546399,
						 3.60555127546399,
						 6.90410505906933,
						 6.90410505906933,
						 6.90410505906933,
						 6.90410505906933,
						 6.90410505906933,
						 6.90410505906933,
						 9.26282894152753,
						 9.26282894152753,
						 9.26282894152753,
						 9.26282894152753,
						 9.26282894152753,
						 9.26282894152753,
						 24.5071418162135,
						 24.5071418162135,
						 24.5071418162135,
						 31.6385840391128,
						 31.6385840391128,
						 31.6385840391128,
						 31.6385840391128,
						 31.6385840391128,
						 31.6385840391128,
						 31.6385840391128,
						 31.6385840391128,
						 31.6385840391128,
						 31.6385840391128,
						 31.6385840391128,
						 31.6385840391128
						}
                    );

        SetCartCNorm(8,
						{1,
						 1,
						 1,
						 3.87298334620742,
						 3.87298334620742,
						 3.87298334620742,
						 3.87298334620742,
						 3.87298334620742,
						 3.87298334620742,
						 8.06225774829855,
						 8.06225774829855,
						 8.06225774829855,
						 8.06225774829855,
						 8.06225774829855,
						 8.06225774829855,
						 11.9582607431014,
						 11.9582607431014,
						 11.9582607431014,
						 11.9582607431014,
						 11.9582607431014,
						 11.9582607431014,
						 13.5593931596198,
						 13.5593931596198,
						 13.5593931596198,
						 35.8747822293042,
						 35.8747822293042,
						 35.8747822293042,
						 35.8747822293042,
						 35.8747822293042,
						 35.8747822293042,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281,
						 54.799635035281
						}
                    );
    }
};



} } // close namespace panache::reorder



#endif
