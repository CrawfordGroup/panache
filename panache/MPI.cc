/*! \file
 *  \brief Functions for handling some mpi stuff (source)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/Parallel.h"
#include "panache/Exception.h"

namespace
{

int * argc_ = nullptr;
char *** argv_ = nullptr;

bool initialized_ = false;

int rank_ = 0;
int size_ = 0;


CTF_World * ctfworld_ = nullptr;

} // close anonymous namespace


namespace panache
{
namespace parallel
{

void Init(int * argc, char *** argv)
{
    if(!initialized_)
    {
        argc_ = argc;
        argv_ = argv;
        initialized_ = true;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);


        #ifdef PANACHE_CYCLOPS
        ctfworld_ = new CTF_World(MPI_COMM_WORLD); 
        #endif
    }
}

void Finalize(void)
{
    if(initialized_)
    {
        #ifdef PANACHE_CYCLOPS
        delete ctfworld_;
        #endif
    }
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
    return size_;    
}

int Rank(void)
{
    return rank_;    
}

bool IsMaster(void)
{
    return (rank_ == 0);
}

std::pair<int, int> MyRange(int totalsize)
{
    int rank = Rank();
    int nelements = (totalsize / Size());
    int start = nelements * rank;

    int leftover = (totalsize % Size());

    if(rank < leftover)
    {
        start += rank;
        nelements++;
    }
    else
        start += leftover;
    
    return std::pair<int, int>(start, start+nelements);
}

#ifdef PANACHE_CYCLOPS
CTF_World & CTFWorld(void)
{
    if(!initialized_)
        throw RuntimeError("MPI is not initialized!");
    return *ctfworld_;
}
#endif

} //close namespace panache
} //close namespace parallel


