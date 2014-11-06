/*! \file
 * \brief Generic, local 4-index tensor storage (header)
 * \ingroup fourindexgroup
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#ifndef PANACHE_LOCAL4INDEXTENSOR_H
#define PANACHE_LOCAL4INDEXTENSOR_H

#include "panache/tensorbackend/FourIndexTensor.h"

#include <vector>

namespace panache
{


class LocalQTensor;  // the component 3-index tensor


class Local4IndexTensor : public FourIndexTensor
{
public:
    Local4IndexTensor(LocalQTensor * left, LocalQTensor * right, size_t nperbatch = 1);

protected:
    virtual int GetNBatches(void) const;
    virtual int GetNLocalIntegrals_(void) const;
    virtual bool GetNextBatch_(void);
    virtual FourIndexIntegral LocalIntegral_(int index);

private:
    LocalQTensor * left_;
    LocalQTensor * right_;

    // for convenience
    int ldim12_;
    int rdim12_;
    int naux_;

    int nperbatch_;
    int nbatches_;
    int curbatch_;

    std::vector<double> contscratch_;
    std::vector<double> integrals_;

    FourIndexTensor::IndexArray DecomposeIndex_(int internalindex) const;

    void GetCurBatch_(void);
};

} // close namespace panache

#endif
