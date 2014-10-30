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
class Molecule;
typedef std::shared_ptr<BasisSet> SharedBasisSet;
typedef std::shared_ptr<Molecule> SharedMolecule;
typedef std::shared_ptr<FittingMetric> SharedFittingMetric;

/*!
 * \brief Generating and manipulation of a density-fitted 3-index tensor
 */
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

    /*!
     * \brief Constructor
     *
     * Initializes the basis set members, and constructs the creates the
     * FittingMetric.
     *
     * \param [in] primary The primary basis set
     * \param [in] auxpath Path to auxiliary basis set file (G98 format)
     * \param [in] directory Full path to a directory to put scratch files
     * \param [in] nthreads Max number of threads to use
     */ 
    DFTensor(SharedBasisSet primary,
             const std::string & auxpath,
             const std::string & directory,
             int nthreads);

protected:
    virtual std::unique_ptr<StoredQTensor> GenQso(int storeflags) const;

private:
    int naux_;   //!< Number of auxiliary basis functions
    SharedBasisSet auxiliary_;  //!< Auxiliary (density fitting) basis set
    SharedFittingMetric fittingmetric_; //!< Fitting metric J

    /// Print the DF tensor information
    void PrintHeader_(void) const;

    /// Initialize some variables
    void Init_(void);

    static SharedBasisSet CreateAuxFromFile_(const std::string & auxpath, SharedMolecule mol);
};

} // close namespace panache

#endif
