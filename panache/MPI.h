/*! \file
 *  \brief Functions for handling some mpi stuff (header)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_MPI_H
#define PANACHE_MPI_H

#include <utility>

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

std::pair<int, int> MyRange(int totalsize);


} //close namespace panache
} //close namespace mpi

#endif

