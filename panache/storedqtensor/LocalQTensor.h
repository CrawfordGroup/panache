/*! \file
 * \brief Generic, local three-index tensor storage (header)
 * \ingroup storedqgroup
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_LOCALQTENSOR_H
#define PANACHE_LOCALQTENSOR_H

#include "panache/storedqtensor/StoredQTensor.h"

namespace panache
{

class TwoBodyAOInt;
typedef std::shared_ptr<TwoBodyAOInt> SharedTwoBodyAOInt;
class FittingMetric;
class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;


/*!
 *  \brief Generic/Abstract interface for locally storing a 3-index tensor
 *  \ingroup storedqgroup
 *
 *  Classes for storing in memory or on disk are derived from this. This class
 *  handles all the transformations, however.
 */
class LocalQTensor : public StoredQTensor
{
public:
    LocalQTensor();

protected:

    /*!
     * \brief Write data with the orbital index as the slowest index
     *
     * \param [in] data Pointer to memory location to put the data
     * \param [in] nij Number of orbital indices to write
     * \param [in] ijstart Starting combined orbital index
     */
    virtual void Write_(double * data, int nij, int ijstart) = 0;


    /*!
     * \brief Write data with the auxiliary index as the slowest index
     *
     * \param [in] data Pointer to memory location to put the data
     * \param [in] nq Number of auxiliary indices to write
     * \param [in] qstart Starting auxiliary index
     */
    virtual void WriteByQ_(double * data, int nq, int qstart) = 0;

    virtual void GenDFQso_(const SharedFittingMetric & fit,
                           const SharedBasisSet primary,
                           const SharedBasisSet auxiliary,
                           int nthreads);

    virtual void GenCHQso_(const SharedBasisSet primary,
                           double delta,
                           int storeflags,
                           int nthreads);

    virtual void Transform_(const std::vector<TransformMat> & left,
                            const std::vector<TransformMat> & right,
                            std::vector<StoredQTensor *> results,
                            int nthreads);

private:
    /*!
     * \brief Compute the cholesky diagonal
     *
     * The size of the \p eris vector is taken to be the number of threads this
     * function can use, which each thread using one TwoBodyAOInt from the vector.
     *
     * \param [in] eris Objects to calculate 4-center integrals
     * \param [in] target Where to put the diagonal information. Should be nso*nso sized
     */
    void ComputeDiagonal_(std::vector<SharedTwoBodyAOInt> & eris, double * target);

    /*!
     * \brief Compute a row for the cholesky Qso
     *
     * The size of the \p eris vector is taken to be the number of threads this
     * function can use, which each thread using one TwoBodyAOInt from the vector.
     *
     * \param [in] eris Objects to calculate 4-center integrals
     * \param [in] row The row to calculate
     * \param [in] target Where to put the row. Should be nso*nso sized
     */
    void ComputeRow_(std::vector<SharedTwoBodyAOInt> & eris, int row, double* target);
};

} // close namespace panache

#endif
