/*! \file
 * \brief Cholesky tensor generation and manipulation (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_CHTENSOR_H
#define PANACHE_CHTENSOR_H

#include "panache/ThreeIndexTensor.h"

namespace panache {


class FittingMetric;
class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;

class CHTensor : public ThreeIndexTensor
{
public:
    /*!
     * \brief Constructor
     *
     * Initializes the basis set members, and constructs the creates the
     * FittingMetric
     *
     * \param [in] primary The primary basis set
     * \param [in] delta Maximum error in the Cholesky procedure
     * \param [in] directory Full path to a directory to put scratch files
     * \param [in] nthreads Max number of threads to use
     */ 
    CHTensor(SharedBasisSet primary,
             double delta,
             const std::string & directory,
             int nthreads);

protected:
    virtual std::unique_ptr<StoredQTensor> GenQso(int storeflags) const;

private:
    double delta_;
};


} // close namespace panache

#endif
