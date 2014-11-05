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


Local4IndexTensor::Local4IndexTensor(LocalQTensor * left, LocalQTensor * right)
                                    : FourIndexTensor(left->ndim1(), left->ndim2(),
                                                      right->ndim1(), right->ndim2()),
                                      left_(left), right_(right)
{
    contscratch_.resize(left->naux()+right->naux());
}

int Local4IndexTensor::GetNLocalIntegrals_(void) const
{
    // only unique!
    return left_->ndim12() * right_->ndim12();
}

FourIndexIntegral Local4IndexTensor::LocalIntegral_(int index)
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

