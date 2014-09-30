#include <fstream>

#include "panache/Exception.h"
#include "panache/ThreeIndexTensor.h"
#include "panache/Parallel.h"
#include "panache/ERI.h"
#include "panache/FittingMetric.h"
#include "panache/BasisSet.h"
    
#ifdef PANACHE_PROFILE
#include "panache/Output.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace panache
{


//////////////////////////////
// CyclopsQTensor
//////////////////////////////
void ThreeIndexTensor::CyclopsQTensor::Reset_(void)
{
    // nothing needed
}

void ThreeIndexTensor::CyclopsQTensor::Read_(double * data, int nij, int ijstart)
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

void ThreeIndexTensor::CyclopsQTensor::ReadByQ_(double * data, int nq, int qstart)
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

void ThreeIndexTensor::CyclopsQTensor::Clear_(void)
{
    tensor_.release();
}

void ThreeIndexTensor::CyclopsQTensor::Init_(void)
{
    //! \todo Symmetry not implemented
    int dimsIJ[3] = {ndim1(), ndim2(), naux()};
    int dimsQ[3] = {naux(), ndim1(), ndim2()};
    int * dims = dimsIJ;

    int symsNS[3] = {NS, NS, NS};
//    int symIJ[3] = {SY, NS, NS};
//    int symQ[3] = {NS, SY, NS};
    int * syms = symsNS;

    if(byq())
        dims = dimsQ;
/*
    if(packed())
    {
        if(byq())
            syms = symQ;
        else
            syms = symIJ;
    }
*/

    tensor_ = std::unique_ptr<CTF_Tensor>(new CTF_Tensor(3, dims, syms, parallel::CTFWorld(), name().c_str()));
}


ThreeIndexTensor::CyclopsQTensor::CyclopsQTensor()
{

}


void ThreeIndexTensor::CyclopsQTensor::DecomposeIndex_(int64_t index, int & i, int & j, int & q)
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

std::unique_ptr<CTF_Matrix>
ThreeIndexTensor::CyclopsQTensor::FillWithMatrix_(double * mat, int nrow, int ncol, int sym, const char * name)
{
    std::unique_ptr<CTF_Matrix> ret(new CTF_Matrix(nrow, ncol, sym, parallel::CTFWorld(), name));

    int64_t np;
    int64_t * idx;
    double * data;

    ret->read_local(&np, &idx, &data);

    for(int64_t n = 0; n < np; n++)
    {
        // note that indicies are in column major order
        // so we have to convert them
        int64_t col = idx[n] / nrow;
        int64_t row = idx[n] - col * nrow;

        data[n] = mat[row * ncol + col];
    }

    ret->write(np, idx, data);

    free(idx);
    free(data);

    return ret;
}


void ThreeIndexTensor::CyclopsQTensor::GenDFQso_(const std::shared_ptr<FittingMetric> & fit,
                                       const SharedBasisSet primary,
                                       const SharedBasisSet auxiliary,
                                       int nthreads)
{
#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
#endif

    // Distributed J matrix
    std::unique_ptr<CTF_Matrix> ctfj(FillWithMatrix_(fit->get_metric(), naux(), naux(), NS, "Jmat"));

    // Now fill up a base matrix
    SharedBasisSet zero(new BasisSet);

    std::vector<std::shared_ptr<TwoBodyAOInt>> eris;
    for(int i = 0; i < nthreads; i++)
        eris.push_back(GetERI(auxiliary, zero, primary, primary));

    int64_t np;
    int64_t *idx;
    double * data;

    // make a tensor of the same size, but don't copy the data
    CTF_Tensor base(*tensor_, false);

    base.read_local(&np, &idx, &data);

    //! \todo better scheduling? We want sequential blocks so
    //        that compute_basisfunction won't miss all the time,
    //        so dynamic is out
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nthreads)
    #endif
    for(int64_t n = 0; n < np; n++)
    {
        int threadnum = 0;
        #ifdef _OPENMP
        threadnum = omp_get_thread_num(); 
        #endif

        int i, j, q;
        DecomposeIndex_(idx[n], i, j, q);
        data[n] = eris[threadnum]->compute_basisfunction(q, 0, i, j);
    }

    #ifdef PANACHE_PROFILE
    unsigned long bmiss = 0;
    for(const auto & it : eris)
        bmiss += it->BMiss();
    output::printf("\n\n***Cyclops GenDFQso : compute_basisfunction misses = %lu\n\n", bmiss);
    #endif

    base.write(np, idx, data);

    free(idx);
    free(data);

    // contract!
    //! \todo check RVO
    if(byq())
        (*tensor_)["iab"] = (*ctfj)["ij"]*base["jab"];
    else
        (*tensor_)["abi"] = (*ctfj)["ij"]*base["abj"];

#ifdef PANACHE_TIMING
    tim.Stop();
    GenTimer().AddTime(tim);
#endif
}


void ThreeIndexTensor::CyclopsQTensor::GenCHQso_(const SharedBasisSet primary,
                                                 double delta,
                                                 int storeflags,
                                                 int nthreads)
{
    throw RuntimeError("NYI");
}


void ThreeIndexTensor::CyclopsQTensor::Transform_(const std::vector<TransformMat> & left,
                                          const std::vector<TransformMat> & right,
                                          std::vector<StoredQTensor *> results,
                                          int nthreads)
{
    // nthreads is ignored

    for(size_t i = 0; i < left.size(); i++)
    {
        #ifdef PANACHE_TIMING
        Timer tim;
        tim.Start();
        #endif

        CyclopsQTensor * rptr = dynamic_cast<CyclopsQTensor *>(results[i]);
        if(rptr == nullptr)
            throw RuntimeError("Cannot store result of Cyclops contraction in non-cyclops object");

        // create CTF Tensors for left and right
        auto ctfleft = FillWithMatrix_(left[i].first, ndim1(), left[i].second, NS, "LEFT");
        auto ctfright = FillWithMatrix_(right[i].first, ndim2(), right[i].second, NS, "RIGHT");

        // we need an intermediate tensor
        int syms[3] = { NS, NS, NS };
        int dimsQ[3] = { naux(), left[i].second, ndim1() };
        int dimsIJ[3] = { left[i].second, ndim1(), naux() };
        int * dims = dimsIJ;
        if(byq())
            dims = dimsQ;

        CTF_Tensor cq(3, dims, syms, parallel::CTFWorld(), "CQ");

        // contract into result tensor
        if(byq())
        {
            if(rptr->byq())
            {
                cq["qij"] = ((*ctfleft)["ai"] * (*tensor_)["qaj"]);
                (*(rptr->tensor_))["qij"] = cq["qia"] * (*ctfright)["aj"];
            }
            else
            {
                cq["qij"] = ((*ctfleft)["ai"] * (*tensor_)["qaj"]);
                (*(rptr->tensor_))["ijq"] = cq["qia"] * (*ctfright)["aj"];
            }
        }
        else
        {
            if(rptr->byq())
            {
                cq["ijq"] = ((*ctfleft)["ai"] * (*tensor_)["ajq"]);
                (*(rptr->tensor_))["qij"] = cq["iaq"] * (*ctfright)["aj"];
            }
            else
            {
                cq["ijq"] = ((*ctfleft)["ai"] * (*tensor_)["ajq"]);
                (*(rptr->tensor_))["ijq"] = cq["iaq"] * (*ctfright)["aj"];
            }
        }

        #ifdef PANACHE_TIMING
        tim.Stop();
        rptr->GenTimer().AddTime(tim);
        #endif
    }
}

} // close namespace panache

