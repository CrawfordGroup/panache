/*! \file
 *  \brief Functions for handling some mpi stuff (header)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_MPI_H
#define PANACHE_MPI_H

#include <utility>

#ifdef PANACHE_CYCLOPS
#include <ctf.hpp>
#endif

namespace panache
{

/*!
 * \brief Holder for some MPI information
 */
namespace mpi
{

void Init(int * argc, char *** argv);

void Finalize(void);

int & Argc(void);
char ** & Argv(void);

int Size(void);
int Rank(void);
bool IsMaster(void);

std::pair<int, int> MyRange(int totalsize);

#ifdef PANACHE_CYCLOPS
CTF_World & CTFWorld(void);
#endif


} //close namespace panache
} //close namespace mpi

#endif

