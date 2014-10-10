/*! \file
 * \brief Three-index tensor storage/manipulation with Cyclops (header)
 * \ingroup storedqgroup
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_CYCLOPSQTENSOR_H
#define PANACHE_CYCLOPSQTENSOR_H

#include <ctf.hpp>

#include "panache/storedqtensor/StoredQTensor.h"

namespace panache
{


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

    /*!
     * \brief Decompose a flattened index into the 3 indicies for a 3-index tensor
     *
     * \param [in] index The flattened (global) index
     * \param [out] i The first orbital index 
     * \param [out] j The second orbital index 
     * \param [out] q The auxiliary index 
     */
    void DecomposeIndex_(int64_t index, int & i, int & j, int & q);

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
    std::unique_ptr<CTF_Matrix> FillWithMatrix_(const double * mat, int nrow, int ncol, int sym, const char * name);


protected:
    virtual void Read_(double * data, int nij, int ijstart);
    virtual void ReadByQ_(double * data, int nq, int qstart);
    virtual void Clear_(void);
    virtual void Init_(void);

    virtual void GenDFQso_(const std::shared_ptr<FittingMetric> & fit,
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

public:
    CyclopsQTensor();
};

} // close namespace panache

#endif


