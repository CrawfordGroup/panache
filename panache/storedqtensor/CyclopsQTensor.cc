/*! \file
 * \brief Three-index tensor storage/manipulation with Cyclops (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/storedqtensor/CyclopsQTensor.h"
#include "panache/Parallel.h"
#include "panache/ERI.h"
#include "panache/FittingMetric.h"
#include "panache/BasisSet.h"
#include "panache/Iterator.h"
#include "panache/Lapack.h"
 
#ifdef PANACHE_PROFILE
#include "panache/Output.h"
#endif

// Using this file requires cyclops & openmp
#include <ctf.hpp>
#include <omp.h>


namespace panache
{


void CyclopsQTensor::Read_(double * data, int nij, int ijstart)
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

void CyclopsQTensor::ReadByQ_(double * data, int nq, int qstart)
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

void CyclopsQTensor::Clear_(void)
{
    tensor_.release();
}

void CyclopsQTensor::Init_(void)
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


CyclopsQTensor::CyclopsQTensor()
{

}


std::unique_ptr<CTF_Matrix>
CyclopsQTensor::FillWithMatrix_(const double * mat, int nrow, int ncol, int sym, const char * name)
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


void CyclopsQTensor::GenDFQso_(const SharedFittingMetric & fit,
                                       const SharedBasisSet primary,
                                       const SharedBasisSet auxiliary,
                                       int nthreads)
{
    int rank = parallel::Rank();
    int maxpershell = primary->max_function_per_shell();
    int maxpershell2 = maxpershell*maxpershell;

    double * J = fit->get_metric();

    // default constructor = zero basis
    SharedBasisSet zero(new BasisSet);

    std::vector<SharedTwoBodyAOInt> eris;
    std::vector<const double *> eribuffers;
    std::vector<double *> A, B;

    int naux = StoredQTensor::naux();

    for(int i = 0; i < nthreads; i++)
    {
        eris.push_back(GetERI(auxiliary, zero, primary, primary));
        eribuffers.push_back(eris.back()->buffer());

        // temporary buffers
        A.push_back(new double[naux*maxpershell2]);
        B.push_back(new double[naux*maxpershell2]);
    }


    const int nauxshell = auxiliary->nshell();
    const int nprimshell = primary->nshell();


    // Get the begin and end shells for this rank
    // Parallelization is over primary basis
    // but ending at shell boundaries
    std::vector<parallel::Range> elranges = parallel::AllRanges((ndim1()*(ndim1()+1))/2);
    std::vector<parallel::Range> shellranges(parallel::Size());
    std::vector<int64_t> allnelements(parallel::Size());

    shellranges[0] = parallel::Range(0,0);

    int64_t nelements = 0;
    int64_t myelements = 0;
    int64_t nsymelements = 0;
    int curproc = 0;
    
    // note: We are calculating just the symmetric part
    // and then applying it to the other triangular part of the
    // tensor as well.
    for(int i = 0, count = 1; i < nprimshell; i++)
    {
        for(int j = 0; j < primary->shell(i).nfunction(); j++, count++)
        {
            nsymelements += count;
            nelements += 2*count-1;
            myelements += 2*count-1;
        }

        if(elranges[curproc].second <= nsymelements)
        {
            shellranges[curproc].second = i+1;
            allnelements[curproc] = myelements;

            curproc++;
            myelements = 0;

            if(curproc < (parallel::Size()))
                shellranges[curproc].first = i+1;
        }
    }

    // allocate. We can reuse the nelements variable
    nelements = allnelements[rank];
    parallel::Range range = shellranges[rank];

    // The above calculation was only on the ij pair.
    // Now we include naux
    nelements *= naux;

    std::unique_ptr<double[]> data(new double[nelements]);
    std::unique_ptr<int64_t[]> idx(new int64_t[nelements]);

    //std::cout << rank << ": MYRANGE: [" << range.first << " , " << range.second << ") NSHELL=" << primary->nshell() << "\n";
    //std::cout << rank << ": NELEMENTS: " << nelements << "\n";

    int64_t curidx = 0;

    for (int M = range.first; M < range.second; M++)
    {
        int threadnum = 0;

#ifdef _OPENMP
        threadnum = omp_get_thread_num();
#endif

        int nm = primary->shell(M).nfunction();
        int mstart = primary->shell(M).function_index();
        int mend = mstart + nm;

        for (int N = 0; N <= M; N++)
        {
            int nn = primary->shell(N).nfunction();
            int nstart = primary->shell(N).function_index();
            int nend = nstart + nn;

            for (int P = 0; P < nauxshell; P++)
            {
                int np = auxiliary->shell(P).nfunction();
                int pstart = auxiliary->shell(P).function_index();
                int pend = pstart + np;

                int ncalc = eris[threadnum]->compute_shell(P,0,M,N);

                if(ncalc)
                {
                    for (int p = pstart, index = 0; p < pend; p++)
                    {
                        for (int m = 0; m < nm; m++)
                        {
                            for (int n = 0; n < nn; n++, index++)
                            {
                                B[threadnum][p*nm*nn + m*nn + n] = eribuffers[threadnum][index];
                            }
                        }
                    }
                }
            }

            // we now have a set of columns of B, although "condensed"
            // we can do a DGEMM with J
            // Access to J are only reads, so that is safe in parallel
            C_DGEMM('T','T',nm*nn, naux, naux, 1.0, B[threadnum], nm*nn, J, naux, 0.0,
                    A[threadnum], naux);


            // write to disk or store in memory
            //! \todo Symmetry exploitation
            if(N == M)
            {
                for (int m0 = 0, m = mstart; m < mend; m0++, m++)
                for (int n0 = 0, n = nstart; n < nend; n0++, n++)
                for (int q = 0; q < naux; q++)
                {
                    idx[curidx] = n*ndim1()*naux + m*naux + q;
                    data[curidx] = A[threadnum][m0*naux*nn + n0*naux + q];
                    curidx++;
                }
            }
            else
            {
                for (int m0 = 0, m = mstart; m < mend; m0++, m++)
                for (int n0 = 0, n = nstart; n < nend; n0++, n++)
                for (int q = 0; q < naux; q++)
                {
                    idx[curidx] = n*ndim1()*naux + m*naux + q;
                    idx[curidx+1] = m*ndim1()*naux + n*naux + q;
                    data[curidx] = data[curidx+1] = A[threadnum][m0*naux*nn + n0*naux + q];
                    curidx += 2;
                }
            }
        }
    }

    //std::cout << rank << " CURIDX/NELEMENTS: " << curidx << " / " << nelements << "\n";

    tensor_->write(curidx, idx.get(), data.get());

    for(int i = 0; i < nthreads; i++)
    {
        delete [] A[i];
        delete [] B[i];
    }

    std::cout << rank << " : Done\n";
}



void CyclopsQTensor::GenCHQso_(const SharedBasisSet primary,
                                                 double delta,
                                                 int storeflags,
                                                 int nthreads)
{
    throw RuntimeError("NYI");
}


void CyclopsQTensor::Transform_(const std::vector<TransformMat> & left,
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

std::pair<int,int>
CyclopsQTensor::ContractMulti_(StoredQTensor * rhs, int ij, int kl, int nij, int nkl,
                               double * out, double * scratch)
{
    CyclopsQTensor * rhsp = dynamic_cast<CyclopsQTensor *>(rhs);
    if(rhsp == nullptr)
        throw RuntimeError("Error - Can't contract CyclopsQTensor with another type!");

    int nq = naux();
    int real_nij = std::min(this->ndim12()-ij, nij);
    int real_nkl = std::min(rhsp->ndim12()-kl, nkl);

/*
    std::array<const int, 3> l_start, l_end;
    std::array<const int, 3> r_start, r_end;

    if(byq())
    {
        ldims = {{0,  }};
    }

    CTF_Tensor lslice = tensor_->slice(

    
    return std::pair<int, int>(real_nij, real_nkl);
*/

    throw RuntimeError("NYI");
}

int CyclopsQTensor::ContractMulti_(StoredQTensor * rhs,
                                   int i, int j, int k, int l, 
                                   int ni, int nj, int nk, int nl, 
                                   double * out, double * scratch)
{

/*
    int nq = naux();
    int nij = ni * nj;
    int nkl = nk * nk;

    int total = nij*nkl;
    if(total == 0)
        return 0; // or throw exception?

    // read (i j | Q)  into scratch
    CyclopsQTensor * rhsp = dynamic_cast<CyclopsQTensor *>(rhs);
    if(rhsp == nullptr)
        throw RuntimeError("Error - Can't contract LocalQTensor with another type!");


    return total;
*/

    throw RuntimeError("NYI");
}


} // close namespace panache

