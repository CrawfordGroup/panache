/*! \file
 * \brief Generic, local 4-index tensor storage (source)
 * \ingroup fourindexgroup
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/tensorbackend/Local4IndexTensor.h"
#include "panache/tensorbackend/LocalQTensor.h"
#include "panache/Exception.h"

namespace panache
{


Local4IndexTensor::Local4IndexTensor(LocalQTensor * left, LocalQTensor * right, size_t nperbatch)
                                    : left_(left), right_(right)
{
    ldim12_ = left->ndim12();
    rdim12_ = right->ndim12();
    naux_ = left->naux();

    if(nperbatch == 0)
        nperbatch_ = ldim12;
    else
        nperbatch_ = std::min(nperbatch, ldim12); // in case nperbatch > the acutal number we have to do

    nbatches_ = ldim12 / nperbatch_;
    if(ldim12 % nperbatch_ > 0)
        nbatches_++;  // last incomplete batch

    contscratch_.resize((rdim12_ + nperbatch_)* naux
    integrals_.resize(nperbatch_ * rdim12_);

    curbatch_ = 0;

    GetCurBatch_();
}

int Local4IndexTensor::GetNBatches(void) const
{
    return nbatches_;
}

int Local4IndexTensor::GetNLocalIntegrals_(void) const
{
    // only unique!
    //return left_->ndim12() * right_->ndim12();
}

bool Local4IndexTensor::GetNextBatch(void)
{
    curbatch_++;
    if(curbatch_ >= nbatches_)
        return false;
    GetCurBatch_();
    return true;
}

FourIndexIntegral Local4IndexTensor::LocalIntegral_(int index)
{
    int internalindex = index + (curbatch_ * right_->ndim12());
}

FourIndexTensor::IndexArray Local4IndexTensor::DecomposeIndex_(int internalindex) const
{
    FourIndexTensor::IndexArray indices;

    indices[0] = internalindex / d123_;
    internalindex -= indices[0] * d123_;

    indices[1] = internalindex / d23_;
    internalindex -= indices[1] * d23_;

    indices[2] = internalindex / dimensions_[3];
    indices[3] = internalindex - (indices[2] * dimensions_[3]);

    return indices;
}

void Local4IndexTensor::GetCurBatch_(void)
{
    // reuse existing infrastructure
    FourIndexTensor::IndexArray ind = DecomposeIndex_(index);

    double val;
    int nint = left_->ContractSingle(right_, ind[0], ind[1], ind[2], ind[3], &val, contscratch_);

    if(nint == 0)
        throw RuntimeError("Error - no integral?");

    int perm = 1;
    if(left_ == right_)
        perm *= 2;

    if(left_->packed() && ind[0] != ind[1])
        perm *= 2;

    if(right_->packed() && ind[2] != ind[3])
        perm *= 2;
        
    return FourIndexIntegral(ind, perm, val);
}

} // close namespace panache

