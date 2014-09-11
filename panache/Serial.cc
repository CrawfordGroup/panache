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

std::pair<int, int> MyRange(int totalsize)
{
    return std::pair<int, int>(0, totalsize);
}

} //close namespace panache
} //close namespace parallel



