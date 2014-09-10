/*! \file
 *  \brief Functions for handling some mpi stuff (source)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/MPI.h"

#include <mpi.h>

namespace
{

int * argc_ = nullptr;
char *** argv_ = nullptr;

bool initialized_ = false;

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
    }
}

void Finalize(void)
{
    if(initialized_)
        MPI_Finalize();
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

} //close namespace panache
} //close namespace mpi


