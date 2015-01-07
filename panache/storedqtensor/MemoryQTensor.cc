/*! \file
 * \brief Three-index tensor storage in memory (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/storedqtensor/MemoryQTensor.h"

namespace panache
{


void MemoryQTensor::Write_(double * data, int nij, int ijstart)
{
    int inaux = naux();

    if(byq())
    {
        int indim12 = ndim12();

        for(int q = 0, qoff = 0; q < inaux; q++, qoff += indim12)
            for(int ij = 0, ijoff = 0; ij < nij; ij++, ijoff += inaux)
                data_[qoff+ijstart+ij] = data[ijoff+q];
    }
    else
    {
        double * start = data_.get() + ijstart * inaux;
        std::copy(data, data+inaux, start);
    }
}


void MemoryQTensor::WriteByQ_(double * data, int nq, int qstart)
{
    int indim12 = ndim12();

    if(byq())
    {
        std::copy(data,
                  data + nq*indim12,
                  data_.get()+qstart*indim12);
    }
    else
    {
        int inaux = naux();

        for(int q = 0, qoff = 0; q < nq; q++, qoff += indim12)
            for(int ij = 0, ijoff = 0; ij < indim12; ij++, ijoff += inaux)
                data_[ijoff+qstart+q] = data[qoff+ij];
    }
}

void MemoryQTensor::Read_(double * data, int nij, int ijstart)
{
    int inaux = naux();

    if(byq())
    {
        int indim12 = ndim12();

        for(int q = 0, qoff = 0; q < inaux; q++, qoff += indim12)
            for(int ij = 0, ijoff = 0; ij < nij; ij++, ijoff += inaux)
                data[ijoff+q] = data_[qoff+ijstart+ij];
    }
    else
    {
        double * start = data_.get() + ijstart * inaux;
        std::copy(start, start+nij*inaux, data);
    }
}

void MemoryQTensor::ReadByQ_(double * data, int nq, int qstart)
{
    int indim12 = ndim12();

    if(byq())
    {
        double * start = data_.get()+qstart*indim12;
        std::copy(start, start+nq*indim12, data);
    }
    else
    {
        int inaux = naux();

        for(int q = 0, qoff = 0; q < nq; q++, qoff += indim12)
            for(int ij = 0, ijoff = 0; ij < indim12; ij++, ijoff += inaux)
                data[qoff+ij] = data_[ijoff+qstart+q];
    }
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

