/*! \file
 *  \brief Functions for managing output to the terminal (header)
 *  \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_OUTPUT_H
#define PANACHE_OUTPUT_H

#include <cstdio>
#include <cstdarg>
#include <ostream>


namespace panache
{
namespace output
{

/*!
 * \brief Sets the output to a C-style FILE pointer
 *
 * \param [in] outfile Where text output should go
 */
void SetOutput(FILE * outfile);


/*!
 * \brief Sets the output to a C++ stream
 *
 * \param [in] outstream Where text output should go
 */
void SetOutput(std::ostream * outstream);


/*!
 * \brief Prints output to that set by SetOutput
 *
 * If output hasn't been set, nothing is output. The formatting
 * is identical to that of the printf family of functions (in fact
 * this function uses those internally).
 *
 * \param [in] format printf-style format string
 * \param [in] ... Arguments for the format string
 */
void printf(const char * format, ...);


} //close namespace panache
} //close namespace output

#endif

