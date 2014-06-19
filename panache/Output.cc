/*! \file
 *  \brief Functions for managing output to the terminal (source)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "Output.h"

namespace
{

#define BUFSIZE 256

FILE * _outfile = NULL;

char _buffer[BUFSIZE];

std::ostream * _outstream = NULL;

bool use_stream = true;

} // close anonymous namespace


namespace panache
{
namespace output
{


void SetOutput(FILE * outfile)
{
    _outfile = outfile;
    use_stream = false;
}

void SetOutput(std::ostream * outstream)
{
    _outstream = outstream;
    use_stream = true;
}


void printf(const char * format, ...)
{
    if(use_stream && _outstream)
    {
        va_list args;
        va_start (args, format);
        vsnprintf(_buffer, BUFSIZE, format, args);
        va_end(args);
        (*_outstream) << _buffer;
    }

    if(!use_stream && _outfile)
    {
        va_list args;
        va_start (args, format);
        vfprintf(_outfile, format, args);
        va_end(args);
    }
}

} //close namespace panache
} //close namespace output


