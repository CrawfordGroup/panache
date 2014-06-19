/*! \file
 * \brief Defines some exception classes 
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_EXCEPTION_H
#define PANACHE_EXCEPTION_H

#include <stdexcept>

namespace panache {

    /*!
     * \brief Panache runtime error exception
     */ 
    typedef std::runtime_error RuntimeError;
}

#endif //PANACHE_EXCEPTION_H
