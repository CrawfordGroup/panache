#include "panache/tensorbackend/FourIndexTensor.h"

namespace panache {

FourIndexTensor::FourIndexTensor(int d1, int d2, int d3, int d4)
                                : dimensions_({{d1, d2, d3, d4}})
{
    // note - arguments 1-based
    // but d123_, etc, refer to zero-based
    d123_ = d2*d3*d4;
    d23_ = d3*d4;

    d012_ = d1*d2*d3;
    d12_ = d2*d3;
}


int FourIndexTensor::GetNLocalIntegrals(void) const
{
    return GetNLocalIntegrals_();
}


FourIndexIntegral FourIndexTensor::LocalIntegral(int index)
{
    return LocalIntegral_(index);
}

FourIndexTensor::IndexArray FourIndexTensor::GetDimensions_(void) const
{
    return dimensions_;
}



FourIndexTensor::IndexArray FourIndexTensor::DecomposeIndex_(int index, bool colmajor) const
{
    FourIndexTensor::IndexArray indices;

    if(colmajor)
    {
        indices[3] = index / d012_;    
        index -= indices[3] * d012_;

        indices[2] = index / d12_;
        index -= indices[2] * d12_;

        indices[1] = index / dimensions_[2];
        indices[0] = index - (indices[1] * dimensions_[2]);
    }
    else
    {
        indices[0] = index / d123_;
        index -= indices[0] * d123_;

        indices[1] = index / d23_;
        index -= indices[1] * d23_;

        indices[2] = index / dimensions_[3];
        indices[3] = index - (indices[2] * dimensions_[3]);
    }

    return indices;
}


} // close namespace panache

