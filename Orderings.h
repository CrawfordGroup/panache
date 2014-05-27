#ifndef PANACHE_REORDER_H
#define PANACHE_REORDER_H

#include <vector>
#include <array>
#include <functional>

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
    array<vector<unsigned short>, MAX_REORDER_AM > orderings_;

public:
    void SetOrder(int am, vector<unsigned short> order)
    {
        if(am > MAX_REORDER_AM)
            throw RuntimeError("Cannot set reordering for this AM - set MAX_REORDER_AM to a larger value");

        orderings_[am].assign(order.begin(), order.end());
    }


    vector<unsigned short> GetOrder(int am) const
    {
        if(am > MAX_REORDER_AM)
            throw RuntimeError("Cannot get reordering for this AM - set MAX_REORDER_AM to a larger value");

        return orderings_[am];
    }


    bool NeedsReordering(int am) const
    {
        if(am > MAX_REORDER_AM)
            throw RuntimeError("Cannot get reordering for this AM - set MAX_REORDER_AM to a larger value");

        return orderings_[am].size();
    }


};




} } // close namespace panache::reorder



#endif
