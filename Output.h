#ifndef PANACHE_OUTPUT_H
#define PANACHE_OUTPUT_H

#include <cstdio>
#include <cstdarg>

#ifndef PANACHE_OUTPUT_C
#include <ostream>
#endif



namespace panache
{
namespace output
{

#ifdef PANACHE_OUTPUT_C
void SetOutput(FILE * outfile);
#else
void SetOutput(std::ostream * outstream);
#endif

void printf(const char * format, ...);

} //close namespace panache
} //close namespace output

#endif

