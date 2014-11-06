#include "panache/tensorbackend/FourIndexTensor.h"

namespace panache {

FourIndexTensor::FourIndexTensor()
{
}

int FourIndexTensor::GetNBatches(void) const
{
    return GetNBatches_();
}

int FourIndexTensor::GetNLocalIntegrals(void) const
{
    return GetNLocalIntegrals_();
}

bool FourIndexTensor::GetNextBatch(void)
{
    return GetNextBatch_();
}

FourIndexIntegral FourIndexTensor::LocalIntegral(int index)
{
    return LocalIntegral_(index);
}


} // close namespace panache

