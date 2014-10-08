/*! \file
 * \brief Generic, local three-index tensor storage (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_LOCALQTENSOR_H
#define PANACHE_LOCALQTENSOR_H

#include "panache/storedqtensor/StoredQTensor.h"

namespace panache
{

class TwoBodyAOInt;
class FittingMetric;
class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet;


/*!
 *  \brief Generic/Abstract interface for locally storing a 3-index tensor
 *
 *  Classes for storing in memory or on disk are derived from this. This class
 *  handles all the transformations, however.
 */
class LocalQTensor : public StoredQTensor
{
public:
    LocalQTensor();

protected:
    virtual void Write_(double * data, int nij, int ijstart) = 0;
    virtual void WriteByQ_(double * data, int nij, int ijstart) = 0;

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

private:
    // for cholesky
    void ComputeDiagonal_(std::vector<std::shared_ptr<TwoBodyAOInt>> & eris, double * target);
    void ComputeRow_(std::vector<std::shared_ptr<TwoBodyAOInt>> & eris, int row, double* target);
};

} // close namespace panache

#endif
