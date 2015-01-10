/*! \file
 * \brief Three-index tensor storage/manipulation with Cyclops (header)
 * \ingroup storedqgroup
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_CYCLOPSQTENSOR_H
#define PANACHE_CYCLOPSQTENSOR_H

#include <ctf.hpp>

#include "panache/storedqtensor/StoredQTensor.h"
#include "panache/Parallel.h"

namespace panache
{


class TwoBodyAOInt;
typedef std::shared_ptr<TwoBodyAOInt> SharedTwoBodyAOInt;
class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet; 


/*!
 *  \brief Class for tensor manipulation using the Cyclops library
 *  \ingroup storedqgroup
 */
class CyclopsQTensor : public StoredQTensor
{
private:
    /// Actual tensor data
    std::unique_ptr<CTF_Tensor> tensor_;

    /// Range of (packed) shell pairs for this process
    parallel::Range myrange_;

    /// Local storage buffer for values
    std::unique_ptr<double[]> mydata_;

    /// Local storage buffer for indicies
    std::unique_ptr<int64_t[]> myidx_;

    /// Number of individual basis function elements in that range
    int64_t mynelements_;

    /*!
     * \brief Convert a locally-stored matrix into a CTF matrix
     *
     * The matrix should be row-major ordered.
     *
     * \note Unsure about AS in \p sym parameter
     *
     * \param [in] mat The matrix to convert
     * \param [in] nrow Number of rows in the matrix
     * \param [in] ncol Number of columns in the matrix
     * \param [in] sym NS for non-symmetric, SY for symmetric, AS for antisymmetric.
     *
     * \param [in] name Some descriptive name
     */
    static std::unique_ptr<CTF_Matrix> FillWithMatrix_(const double * mat, int nrow, int ncol, int sym, const char * name);


    /*!
     * \brief Find the maximum value and index of a distributed cyclops vector
     *
     * I can't just do vec.reduce(CTF_OP_MAX) since I also need the index (for 
     * Cholesky decomposition).
     *
     * \param [in] vec Vector to find the maximum value in
     * \return The maximum value and corresponding index
     */
    static std::pair<int64_t, double> FindVecMax_(CTF_Vector & vec);


    /*!
     * \brief Compute the cholesky diagonal
     *
     * The size of the \p eris vector is taken to be the number of threads this
     * function can use, which each thread using one TwoBodyAOInt from the vector.
     *
     * \note Not static like in LocalQTensor since we need access to range information
     *       (mynelements_ and myrange_)
     *
     * \param [in] eris Objects to calculate 4-center integrals
     * \param [in] target Where to put the diagonal information. Should be nso*nso sized
     */
    void ComputeDiagonal_(std::vector<SharedTwoBodyAOInt> & eris, CTF_Vector & target);

    /*!
     * \brief Compute a row for the cholesky Qso
     *
     * The size of the \p eris vector is taken to be the number of threads this
     * function can use, which each thread using one TwoBodyAOInt from the vector.
     *
     * \note Not static like in LocalQTensor since we need access to range information
     *       (mynelements_ and myrange_)
     *
     * \param [in] eris Objects to calculate 4-center integrals
     * \param [in] row The row to calculate
     * \param [in] target Where to put the row. Should be nso*nso sized
     */
    void ComputeRow_(std::vector<SharedTwoBodyAOInt> & eris, int64_t row, CTF_Vector & target);

    /*!
     * \brief Calculates the range of shells to calculate for a specific
     *        process.
     *
     * This function aims to calculate as even a distribution as possible
     * for a calculation over basis function pairs with packed indices.
     */
    static std::pair<int64_t, parallel::Range>
    ShellRange2_(const SharedBasisSet & basis);

protected:
    virtual void Read_(double * data, int nij, int ijstart);
    virtual void ReadByQ_(double * data, int nq, int qstart);
    virtual void Init_(void);

    virtual void GenDFQso_(const SharedFittingMetric & fit,
                           const SharedBasisSet primary,
                           const SharedBasisSet auxiliary,
                           int nthreads);

    virtual void GenCHQso_(const SharedBasisSet primary,
                           double delta,
                           int nthreads);

    virtual void Transform_(const std::vector<TransformMat> & left,
                            const std::vector<TransformMat> & right,
                            std::vector<StoredQTensor *> results,
                            int nthreads);

public:
    /*
     * \brief Construct with some basic information
     *
     * \param [in] storeflags How the tensor should be stored (packed, etc)
     * \param [in] name Some descriptive name
     */
    CyclopsQTensor(int storeflags, const std::string & name);
};

} // close namespace panache

#endif


