#include "panache/storedqtensor/MemoryQTensor.h"

namespace panache
{


void MemoryQTensor::Reset_(void)
{
    // nothing needed
}


void MemoryQTensor::Write_(double * data, int nij, int ijstart)
{
    if(byq())
    {
        for(int q = 0; q < naux(); q++)
        for(int ij = 0; ij < nij; ij++)
            data_[q*ndim12()+ijstart+ij] = data[ij*naux()+q];
    }
    else
    {
        double * start = data_.get() + ijstart * naux();
        std::copy(data, data+naux(), start);
    }
}


void MemoryQTensor::WriteByQ_(double * data, int nq, int qstart)
{
    if(byq())
    {
        std::copy(data,
                  data + nq*ndim12(),
                  data_.get()+qstart*ndim12());
    }
    else
    {
        for(int q = 0; q < nq; q++)
        {
            for(int ij = 0; ij < ndim12(); ij++)
                data_[ij*naux()+qstart+q] = data[q*ndim12()+ij];
        }
    }
}

void MemoryQTensor::Read_(double * data, int nij, int ijstart)
{
    // index ij is given by calling function and takes into account packing
    if(byq())
    {
        for(int q = 0; q < naux(); q++)
        for(int ij = 0; ij < nij; ij++)
            data[ij*naux()+q] = data_[q*ndim12()+ijstart+ij];
    }
    else
    {
        double * start = data_.get() + ijstart * naux();
        std::copy(start, start+nij*naux(), data);
    }
}

void MemoryQTensor::ReadByQ_(double * data, int nq, int qstart)
{
    if(byq())
    {
        double * start = data_.get()+qstart*ndim12();
        std::copy(start, start+nq*ndim12(), data);
    }
    else
    {
        for(int q0 = 0, q = qstart; q0 < nq; q++, q0++)
        for(int m = 0; m < ndim1(); m++)
        for(int n = 0; n < ndim2(); n++)
            data[q0*ndim12() + calcindex(m,n)]
                = data_[calcindex(m,n)*naux() + q];
    }
}

void MemoryQTensor::Clear_(void)
{
    data_.reset();
}

void MemoryQTensor::Init_(void)
{
    if(!data_)
        data_ = std::unique_ptr<double []>(new double[storesize()]);
}

MemoryQTensor::MemoryQTensor()
{
}


} // close namespace panache

