#include <fstream>

#include "panache/Exception.h"
#include "panache/DFTensor.h"
#include "panache/Parallel.h"
#include "panache/ERI.h"
#include "panache/FittingMetric.h"
#include "panache/BasisSet.h"

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

    tensor_.read(nelements, indices.data(), data);
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

    tensor_.read(nelements, indices.data(), data);
}

void DFTensor::CyclopsQTensor::Clear_(void)
{
    tensor_ = CTF_Tensor();
}

void DFTensor::CyclopsQTensor::Init_(void)
{
}


DFTensor::CyclopsQTensor::CyclopsQTensor(int naux, int ndim1, int ndim2, int storeflags, const std::string & name)
            : StoredQTensor(naux, ndim1, ndim2, storeflags), name_(name)
{
}


void DFTensor::CyclopsQTensor::DecomposeIndex_(int64_t index, int & i, int & j, int & q)
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


CTF_Matrix DFTensor::CyclopsQTensor::FillWithMatrix_(double * mat, int nrow, int ncol, int sym, const char * name)
{
    CTF_Matrix ret(nrow, ncol, sym, parallel::CTFWorld(), name);

    int64_t np;
    int64_t * idx;
    double * data;

    ret.read_local(&np, &idx, &data);

    for(int64_t n = 0; n < np; n++)
    {
        // note that indicies are in column major order
        // so we have to convert them
        int64_t col = idx[n] / nrow;
        int64_t row = idx[n] - col * nrow;

        data[n] = mat[row * ncol + col];
    }

    ret.write(np, idx, data);

    free(idx);
    free(data);

    return ret;
}


void DFTensor::CyclopsQTensor::GenQso_(const std::shared_ptr<FittingMetric> & fit,
                                       const SharedBasisSet primary,
                                       const SharedBasisSet auxiliary,
                                       int nthreads)
{
    // nthreads is ignored

    int dimsIJ[3] = {ndim1(), ndim2(), naux()};
    int dimsQ[3] = {naux(), ndim1(), ndim2()};
    int * dims = dimsIJ;

    int symsNS[3] = {NS, NS, NS};
    int symIJ[3] = {SY, NS, SY};
    int symQ[3] = {NS, SY, NS};
    int * syms = symsNS;

    if(byq())
        dims = dimsQ;

    if(packed())
    {
        if(byq())
            syms = symQ;
        else
            syms = symIJ;
    }


    // Distributed J matrix
    CTF_Matrix ctfj = FillWithMatrix_(fit->get_metric(), naux(), naux(), NS, "Jmat");

    // Now fill up a base matrix
    SharedBasisSet zero(new BasisSet);
    std::shared_ptr<TwoBodyAOInt> eri = GetERI(auxiliary, zero, primary, primary);

    int64_t np;
    int64_t *idx;
    double * data;

    CTF_Tensor base(3, dims, syms, parallel::CTFWorld(), name_.c_str());

    base.read_local(&np, &idx, &data);

    int i, j, q;

    //! \todo openmp here?
    for(int64_t n = 0; n < np; n++)
    {
        DecomposeIndex_(idx[n], i, j, q);
        data[n] = eri->compute_basisfunction(q, 0, i, j);
    }

    base.write(np, idx, data);

    free(idx);
    free(data);

    // contract!
    if(byq())
        tensor_["iab"] = ctfj["ij"]*base["jab"];
    else
        tensor_["abi"] = ctfj["ij"]*base["abj"];
}


void DFTensor::CyclopsQTensor::Transform_(const std::vector<TransformMat> & left,
                                          const std::vector<TransformMat> & right,
                                          std::vector<StoredQTensor *> results,
                                          int nthreads)
{
    // nthreads is ignored

    // fi
    
}


} // close namespace panache

