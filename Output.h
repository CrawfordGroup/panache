#ifndef PANACHE_OUTPUT_H
#define PANACHE_OUTPUT_H

#include <cstdio>
#include <cstdarg>
#include <ostream>


namespace panache
{
namespace output
{

void SetOutput(FILE * outfile);
void SetOutput(std::ostream * outstream);

void printf(const char * format, ...);

} //close namespace panache
} //close namespace output

#endif

