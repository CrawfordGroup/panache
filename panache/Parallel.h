/*! \file
 *  \brief Functions for handling some mpi stuff (header)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_PARALLEL_H
#define PANACHE_PARALLEL_H

#include <utility>
#include <vector>
#include <cstdint>

#ifdef PANACHE_CYCLOPS
#include <ctf.hpp>
#endif

namespace panache
{

/*!
 * \brief Holder for some MPI information
 */
namespace parallel
{

void Init(int * argc, char *** argv);

void Finalize(void);

int & Argc(void);
char ** & Argv(void);

int Size(void);
int Rank(void);
bool IsMaster(void);

typedef std::pair<int64_t, int64_t> Range;

Range MyRange(int64_t totalsize);
std::vector<Range> AllRanges(int64_t totalsize);

#ifdef PANACHE_CYCLOPS
CTF_World & CTFWorld(void);
#endif


} //close namespace panache
} //close namespace parallel

#endif

