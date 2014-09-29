/*! \file
 * \brief Density fitting tensor generation and manipulation (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_DFTENSOR_H
#define PANACHE_DFTENSOR_H

#include "panache/ThreeIndexTensor.h"

namespace panache {

class FittingMetric;
class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;


class DFTensor : public ThreeIndexTensor
{
public:
    /*!
     * \brief Constructor
     *
     * Initializes the basis set members, and constructs the creates the
     * FittingMetric
     *
     * \param [in] primary The primary basis set
     * \param [in] auxiliary The auxiliary (DF) basis set
     * \param [in] directory Full path to a directory to put scratch files
     * \param [in] nthreads Max number of threads to use
     */ 
    DFTensor(SharedBasisSet primary,
             SharedBasisSet auxiliary,
             const std::string & directory,
             int nthreads);

protected:
    virtual std::unique_ptr<StoredQTensor> GenQso(int storeflags) const;

private:
    int naux_;   //!< Number of auxiliary basis functions
    SharedBasisSet auxiliary_;
    std::shared_ptr<FittingMetric> fittingmetric_; //!< Fitting metric J
};

} // close namespace panache

#endif
