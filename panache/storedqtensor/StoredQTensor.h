/*! \file
 * \brief Generic three-index tensor storage (header)
 * \ingroup storedqgroup
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef PANACHE_STOREDQTENSOR_H
#define PANACHE_STOREDQTENSOR_H

#include <memory>
#include <string>
#include <vector>

#include "panache/Timing.h"

namespace panache
{

class FittingMetric;
typedef std::shared_ptr<FittingMetric> SharedFittingMetric;
class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;


/*!
 *  \brief Generic/Abstract interface for storing a 3-index tensor
 *  \ingroup storedqgroup
 */
class StoredQTensor
{
public:

    /*
     * \brief Construct with some basic information
     *
     * \param [in] storeflags How the tensor should be stored (packed, etc)
     * \param [in] name Some descriptive name
     */
    StoredQTensor(int storeflags, const std::string & name);

    virtual ~StoredQTensor();

    /*! 
     * \brief Get what the current storage flags are
     * 
     * See Flags.h
     */
    int StoreFlags(void) const;

   
    /*!
     * \brief Read data with the orbital index as the slowest index
     *
     * The buffer must be large enough to hold the requested data.
     * This is the caller's responsibility - no checking is done at this point.
     * The buffer should be nij * naux elements.
     *
     * \param [in] data Pointer to memory location to put the data
     * \param [in] nij Number of orbital pairs to obtain
     * \param [in] ijstart Starting orbital index
     */
    int Read(double * data, int nij, int ijstart);

   
    /*!
     * \brief Read data with the auxiliary index as the slowest index
     *
     * The buffer must be large enough to hold the requested data.
     * This is the caller's responsibility - no checking is done at this point.
     * The buffer should be nq * ndim12 elements. The returned data will be packed
     * if the tensor is stored that way.
     *
     * \param [in] data Pointer to memory location to put the data
     * \param [in] nq Number of auxiliary indices to obtain
     * \param [in] qstart Starting auxiliary index
     */
    int ReadByQ(double * data, int nq, int qstart);


    /*!
     * \brief Initialize storage for a given size
     *
     * \param [in] naux Size along the auxiliary index
     * \param [in] ndim1 Size along the first dimension
     * \param [in] ndim2 Size along the second dimension
     */
    void Init(int naux, int ndim1, int ndim2);


    /*!
     * \brief Initialize storage given the sizes of another StoredQTensor
     */
    void Init(const StoredQTensor & rhs);
  

    /*!
     * \brief Generate density-fitted Qso tensor
     *
     * \param [in] fit Pre-calculated fitting metric
     * \param [in] primary Primary basis set
     * \param [in] auxiliary Auxiliary basis set
     * \param [in] nthreads Number of threads to use
     */ 
    void GenDFQso(const SharedFittingMetric & fit,
                  const SharedBasisSet primary,
                  const SharedBasisSet auxiliary,
                  int nthreads);

    /*!
     * \brief Generate cholesky-based Qso tensor
     *
     * \param [in] primary Primary basis set
     * \param [in] delta Maximum error in the cholesky procedure
     * \param [in] nthreads Number of threads to use
     */ 
    void GenCHQso(const SharedBasisSet primary,
                  double delta,
                  int nthreads);

    /*!
     * \brief A matrix used to transform Qso.
     *
     * The first member is a pointer to the matrix, and the
     * second member is the size along the second dimension.
     * That is, the matrix is expected to be nso x TransformMat.second.
     */
    typedef std::pair<double *, int> TransformMat;


    /*!
     *  \brief Transform this tensor
     *
     *  This function carries out multiple transformations. First,
     *  it reads by q, and then for each TransformMat in \p left and \p right,
     *  stores it in the StoredQTensor object in \p results.
     *
     *  The objects in the \p results vector must be created and initialized already.
     */
    void Transform(const std::vector<TransformMat> & left,
                   const std::vector<TransformMat> & right,
                   std::vector<StoredQTensor *> results,
                   int nthreads);


    /// Get the timer for generation of this tensor
    CumulativeTime & GenTimer(void);

    /// Get the timer for getting batches by IJ
    CumulativeTime & GetBatchTimer(void);

    /// Get the timer for getting batches by Q
    CumulativeTime & GetQBatchTimer(void);


    /// Get the size along the auxiliary index
    int naux(void) const;

    /// Get the size along the first dimension
    int ndim1(void) const;

    /// Get the size along the second dimension
    int ndim2(void) const;

    /// Get the size with respect to combined orbital indices. Depends on packing.
    int ndim12(void) const;

    /// Get all storage flags
    int storeflags(void) const;

    /// Get whether or not this tensor is stored packed
    int packed(void) const;

    /// Get whether or not this tensor is stored by q
    int byq(void) const;

    /// Get whether or not this tensor is filled in
    bool filled(void) const;

    /*!
     *  \brief Calculate the combined orbital index given a pair of orbitals.
     *
     * Depends on packing. If packed, returns i*(i+1)/2+j. Else, returns i*ndim2 + j.
     */
    int calcindex(int i, int j) const;

    /// Get the name of this tensor
    const std::string & name(void) const;

protected:

    /*!
     * \brief Initialize the object
     *
     * To be implemented by derived classes
     *
     * When this is called, all the dimensions in the base class have
     * been filled in
     */
    virtual void Init_(void) = 0;

    /// \copydoc Read()
    /// To be implemented by derived classes
    virtual void Read_(double * data, int nij, int ijstart) = 0;

    /// \copydoc ReadByQ()
    /// To be implemented by derived classes
    virtual void ReadByQ_(double * data, int nq, int qstart) = 0;

    /// \copydoc GenDFQso()
    /// To be implemented by derived classes
    virtual void GenDFQso_(const SharedFittingMetric & fit,
                           const SharedBasisSet primary,
                           const SharedBasisSet auxiliary,
                           int nthreads) = 0;

    /// \copydoc GenCHQso()
    /// To be implemented by derived classes
    virtual void GenCHQso_(const SharedBasisSet primary,
                           double delta,
                           int nthreads) = 0;


    /// \copydoc Transform()
    /// To be implemented by derived classes
    virtual void Transform_(const std::vector<TransformMat> & left,
                            const std::vector<TransformMat> & right,
                            std::vector<StoredQTensor *> results,
                            int nthreads) = 0;

    /// Get the total size of the stored tensor
    int storesize(void) const;

    /// Mark this tensor object as filled in
    void markfilled(void);

private:
    std::string name_;

    int naux_;   //!< Number of auxiliary functions
    int ndim1_;  //!< Length of index 1
    int ndim2_;  //!< Length of index 2
    int ndim12_; //!< Combined size of index 1 and 2 (depends on packing)
    int storeflags_; //!< How is this tensor stored
    bool filled_;    //!< Tensor has been filled in

    CumulativeTime gen_timer_; //!< Timer for the generation of this tensor
    CumulativeTime getijbatch_timer_; //!< Timer for getting batch by orbital index
    CumulativeTime getqbatch_timer_; //!< Timer for getting batch by q index

    // disable copying, etc
    StoredQTensor & operator=(const StoredQTensor & rhs) = delete;
    StoredQTensor(const StoredQTensor & rhs) = delete;
    StoredQTensor(const StoredQTensor && rhs) = delete;
};

typedef std::unique_ptr<StoredQTensor> UniqueStoredQTensor;

} // close namespace panache

#endif

