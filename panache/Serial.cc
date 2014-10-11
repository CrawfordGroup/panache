/*! \file
 *  \brief Serial parallelization fallback (source)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/Parallel.h"

namespace
{

int * argc_ = nullptr;
char *** argv_ = nullptr;

}

namespace panache
{
namespace parallel
{

void Init(int * argc, char *** argv)
{
}

void Finalize(void)
{
}

int & Argc(void)
{
    return *argc_;
}

char ** & Argv(void)
{
    return *argv_;
}

int Size(void)
{
    return 1;
}

int Rank(void)
{
    return 0;
}

bool IsMaster(void)
{
    return true;
}

Range MyRange(int64_t totalsize)
{
    return Range(0, totalsize);
}

std::vector<Range> AllRanges(int64_t totalsize)
{
    std::vector<Range> ret;
    ret.push_back(Range(0, totalsize));
    return ret;
}

} //close namespace panache
} //close namespace parallel



