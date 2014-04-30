#include "Output.h"

namespace
{

#ifdef PANACHE_OUTPUT_C
FILE * _outfile = NULL;
#else
#define BUFSIZE 256
char _buffer[BUFSIZE];
std::ostream * _outstream = NULL;
#endif

} // close anonymous namespace


namespace panache
{
namespace output
{

#ifdef PANACHE_OUTPUT_C

void SetOutput(FILE * outfile)
{
    _outfile = outfile;
}

void printf(const char * format, ...)
{
    if(_outfile)
    {
        va_list args;
        va_start (args, format);
        vfprintf(_outfile, format, args);
        va_end(args);
    }
}

#else

void SetOutput(std::ostream * outstream)
{
    _outstream = outstream;
}

void printf(const char * format, ...)
{
    if(_outstream)
    {
        va_list args;
        va_start (args, format);
        vsnprintf(_buffer, BUFSIZE, format, args);
        va_end(args);
        (*_outstream) << _buffer;
    }
}

#endif


} //close namespace panache
} //close namespace output

