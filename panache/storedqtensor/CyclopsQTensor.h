/*! \file
 * \brief Three-index tensor storage in memory (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_CYCLOPSQTENSOR_H
#define PANACHE_CYCLOPSQTENSOR_H

#include <memory>
#include <vector>

#include <ctf.hpp>

#include "panache/storedqtensor/StoredQTensor.h"

namespace panache
{


class BasisSet;
typedef std::shared_ptr<BasisSet> SharedBasisSet; 


/*!
 *  \brief Class for itensor manipulation using the Cyclops library
 */
class CyclopsQTensor : public StoredQTensor
{
private:
    std::unique_ptr<CTF_Tensor> tensor_;

    void DecomposeIndex_(int64_t index, int & i, int & j, int & q);
    std::unique_ptr<CTF_Matrix> FillWithMatrix_(double * mat, int nrow, int ncol, int sym, const char * name);


protected:
    virtual void Reset_(void);
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


