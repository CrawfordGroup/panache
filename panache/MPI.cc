/*! \file
 *  \brief Functions for handling some mpi stuff (source)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/MPI.h"
#include "panache/Exception.h"

namespace
{

int * argc_ = nullptr;
char *** argv_ = nullptr;

bool initialized_ = false;

CTF_World * ctfworld_ = nullptr;

} // close anonymous namespace


namespace panache
{
namespace mpi
{

void Init(int * argc, char *** argv)
{
    if(!initialized_)
    {
        argc_ = argc;
        argv_ = argv;
        MPI_Init(argc, argv);
        initialized_ = true;

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

        MPI_Finalize();
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
    int size = 0;

    if(initialized_)
        MPI_Comm_size(MPI_COMM_WORLD, &size);

    return size;    
}

int Rank(void)
{
    int rank = 0;

    if(initialized_)
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    return rank;    
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
} //close namespace mpi


