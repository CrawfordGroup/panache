#include <fstream>

#include "panache/Exception.h"
#include "panache/DFTensor.h"
#include "panache/MPI.h"

namespace panache
{


//////////////////////////////
// CyclopsQTensor
//////////////////////////////
void DFTensor::CyclopsQTensor::Reset_(void)
{
    // nothing needed
}

void DFTensor::CyclopsQTensor::Read_(double * data, int nij, int ijstart)
{
    size_t nelements = nij*naux();

    std::vector<long> indices;
    indices.reserve(nelements);


    IJIterator it(ndim1(), ndim2(), packed());
    it += ijstart;

    if(byq())
    {
        for(int ij = 0; ij < nij; ij++)
        {
            for(int q = 0; q < naux(); q++)
                indices.push_back(q+it.i()*naux()+it.j()*naux()*ndim1());
            ++it;
        }
    }
    else
    {
        for(int ij = 0; ij < nij; ij++)
        {
            for(int q = 0; q < naux(); q++)
                indices.push_back(it.i()+it.j()*ndim1()+q*ndim1()*ndim2());
            ++it;
        }
    }

    tensor_->read(nelements, indices.data(), data);
}

void DFTensor::CyclopsQTensor::ReadByQ_(double * data, int nq, int qstart)
{
    size_t nelements = nq*ndim12();

    std::vector<long> indices;
    indices.reserve(nelements);


    if(byq())
    {
        for(int q = 0; q < naux(); q++)
        {
            IJIterator it(ndim1(), ndim2(), packed());
            
            for(int ij = 0; ij < ndim12(); ij++)
            {
                indices.push_back((q+qstart)+it.i()*naux()+it.j()*naux()*ndim1());
                ++it;
            }
        }
    }
    else
    {
        for(int q = 0; q < naux(); q++)
        {
            IJIterator it(ndim1(), ndim2(), packed());
            
            for(int ij = 0; ij < ndim12(); ij++)
            {
                indices.push_back(it.i()+it.j()*ndim1()+(q+qstart)*ndim1()*ndim2());
                ++it;
            }
        }
    }

    tensor_->read(nelements, indices.data(), data);
}

void DFTensor::CyclopsQTensor::Clear_(void)
{
    tensor_.release();
}

void DFTensor::CyclopsQTensor::Init_(void)
{
}


DFTensor::CyclopsQTensor::CyclopsQTensor(int naux, int ndim1, int ndim2, int storeflags, const std::string & name)
            : StoredQTensor(naux, ndim1, ndim2, storeflags)
{

    int dims1[3] = {ndim1, ndim2, naux};
    int dims2[3] = {naux, ndim1, ndim2};
    int * dims = dims1;

    int syms1[3] = {NS, NS, NS};
    int syms2[3] = {NS, SY, SY};
    int * syms = syms1;

    if(byq())
        dims = dims2;

    if(packed())
        syms = syms2;

    tensor_ = std::unique_ptr<CTF_Tensor>(new CTF_Tensor(3, dims, syms, mpi::CTFWorld(), name.c_str()));
}


void DFTensor::CyclopsQTensor::DecomposeIndex_(int index, int & i, int & j, int & q)
{
    if(byq())
    {
        j = index / (ndim1() * naux());
        index -= j * ndim1() * naux();
        i = index / naux();
        q = index - i * naux();
    }
    else
    {
        q = index / (ndim1() * ndim2());
        index -= q * ndim1() * ndim2();
        j = index / ndim1();
        i = index - j * ndim1();
    }
}

void DFTensor::CyclopsQTensor::GenQso_(const std::shared_ptr<FittingMetric> & fit,
                                       const SharedBasisSet primary,
                                       const SharedBasisSet auxiliary,
                                       int nthreads)
{
    // nthreads is ignored
}


void DFTensor::CyclopsQTensor::Transform_(const std::vector<TransformMat> & left,
                                          const std::vector<TransformMat> & right,
                                          std::vector<StoredQTensor *> results,
                                          int nthreads)
{
    // nthreads is ignored
}


} // close namespace panache

