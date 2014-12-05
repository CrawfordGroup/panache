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
    //int rank = parallel::Rank();
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


    // Get the begin and end shells for this rank
    // Parallelization is over primary basis
    // but ending at shell boundaries
    auto shellrangeinfo = ShellRange2_(primary);

    int64_t nelements = shellrangeinfo.first;
    parallel::Range range = shellrangeinfo.second;

    // The above calculation was only on the ij pair.
    // Now we include naux
    nelements *= naux;

    std::unique_ptr<double[]> data(new double[nelements]);
    std::unique_ptr<int64_t[]> idx(new int64_t[nelements]);

    //std::cout << rank << ": MYRANGE: [" << range.first << " , " << range.second << ") NSHELL=" << primary->nshell() << "\n";
    //std::cout << rank << ": NELEMENTS: " << nelements << "\n";

    int64_t curidx = 0;

//! \todo Fix threading in MPI?
//#ifdef _OPENMP
//    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//#endif
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

    //std::cout << rank << " : Done\n";
}

std::pair<int64_t, parallel::Range>
CyclopsQTensor::ShellRange2_(const SharedBasisSet & basis)
{
    int nbf = basis->nbf();
    int nbf12 = (nbf*(nbf+1))/2;

    // Get the begin and end shells for this rank
    // Parallelization is over basis
    // but ending at shell boundaries
    std::vector<parallel::Range> elranges = parallel::AllRanges(nbf12);
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
    const int nshell = basis->nshell();
    for(int i = 0, count = 1; i < nshell; i++)
    {
        for(int j = 0; j < basis->shell(i).nfunction(); j++, count++)
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

    int rank = parallel::Rank();
    return std::pair<int64_t, parallel::Range>(allnelements[rank], shellranges[rank]);
}

void CyclopsQTensor::ComputeDiagonal_(std::vector<SharedTwoBodyAOInt> & eris, 
                                      CTF_Vector & target)
{
    SharedBasisSet basis = eris[0]->basis();

    const int nbf = basis->nbf();

    std::unique_ptr<double[]> data(new double[mynelements_]);
    std::unique_ptr<int64_t[]> idx(new int64_t[mynelements_]);

    int64_t curidx = 0;

    //size_t nthreads = eris.size();
//! \todo Fix threading in MPI?
//#ifdef _OPENMP
//    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//#endif
    for (int M = myrange_.first; M < myrange_.second; M++) 
    {
        int threadnum = 0;

#ifdef _OPENMP
        threadnum = omp_get_thread_num();
#endif
        TwoBodyAOInt * integral = eris[threadnum].get();
        const double* buffer = integral->buffer();

        int nM = basis->shell(M).nfunction();
        int mstart = basis->shell(M).function_index();

        for (int N = 0; N <= M; N++) 
        {
            int nint = integral->compute_shell(M,N,M,N);

            if(nint)
            {
                int nN = basis->shell(N).nfunction();
                int nstart = basis->shell(N).function_index();
    
                if (N == M)
                {
                  for (int om = 0; om < nM; om++)
                  {
                    for (int on = 0; on < nN; on++)
                    {
                        idx[curidx] = (om + mstart) * nbf + (on + nstart);
                        data[curidx] = buffer[om * nN * nM * nN + on * nM * nN + om * nN + on];
                        curidx += 1;
                    }
                  }
                }
                else
                {    
                  for (int om = 0; om < nM; om++)
                  {
                    for (int on = 0; on < nN; on++)
                    {
                        idx[curidx] = (om + mstart) * nbf + (on + nstart);
                        data[curidx] = buffer[om * nN * nM * nN + on * nM * nN + om * nN + on];
                        idx[curidx+1] = (on + nstart) * nbf + (om + mstart); 
                        data[curidx+1] = data[curidx];
                        curidx += 2;
                    }
                  }
                }
            }
        }
    }

    target.write(curidx, idx.get(), data.get());
}

void CyclopsQTensor::ComputeRow_(std::vector<SharedTwoBodyAOInt> & eris, 
                                 int row, CTF_Vector & target)
{
    SharedBasisSet basis = eris[0]->basis();

    const int nbf = basis->nbf();

    int r = row / nbf;
    int s = row % nbf;
    int R = basis->function_to_shell(r);
    int S = basis->function_to_shell(s);

    int nR = basis->shell(R).nfunction();
    int nS = basis->shell(S).nfunction();
    int rstart = basis->shell(R).function_index();
    int sstart = basis->shell(S).function_index();

    int oR = r - rstart;
    int os = s - sstart;

/*
    
    //size_t nthreads = eris.size();
//! \todo Fix threading in MPI?
//#ifdef _OPENMP
//    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//#endif
    for (int M = myrange_.first; M < myrange_.second; M++)
    {
        int threadnum = 0;

#ifdef _OPENMP
        threadnum = omp_get_thread_num();
#endif

        TwoBodyAOInt * integral = eris[threadnum].get();
        const double* buffer = integral->buffer();

        int nM = basis->shell(M).nfunction();
        int mstart = basis->shell(M).function_index();

        for (int N = 0; N <= M; N++)
        {
            int nint = integral->compute_shell(M,N,R,S);

            if(nint)
            {
                int nN = basis->shell(N).nfunction();
                int nstart = basis->shell(N).function_index();
    
                for (int om = 0; om < nM; om++) {
                    for (int on = 0; on < nN; on++) {
                        target[(om + mstart) * nbf + (on + nstart)] =
                        target[(on + nstart) * nbf + (om + mstart)] =
                            buffer[om * nN * nR * nS + on * nR * nS + oR * nS + os];
                    }
                }
            }
        }
    }
*/
}

void CyclopsQTensor::GenCHQso_(const SharedBasisSet primary,
                                                 double delta,
                                                 int storeflags,
                                                 int nthreads)
{
    auto shellrangeinfo = ShellRange2_(primary);

    mynelements_ = shellrangeinfo.first;
    myrange_ = shellrangeinfo.second;

    std::vector<SharedTwoBodyAOInt> eris;

    for(int i = 0; i < nthreads; i++)
        eris.push_back(GetERI(primary, primary, primary, primary));

    int nQ = 0;
    int n = primary->nbf();
    int n2 = n*n;


    CTF_Vector diag(n2, parallel::CTFWorld(), "CHODIAG");
    ComputeDiagonal_(eris, diag);

    auto max = FindVecMax_(diag);

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

std::pair<int64_t, double> CyclopsQTensor::FindVecMax_(CTF_Vector & vec)
{
  int myRank = parallel::Rank();
  int numPes = parallel::Size();

  // find maximum on this process
  int64_t np;
  int64_t * idx;
  double * data;

  // \todo I wish this didn't create copies...
  // Or I could just get a list of indicies
  vec.read_local(&np, &idx, &data);

  double max = data[0];
  int64_t index = 0;

  for(int64_t i = 1; i < np; i++)
  {
    if(max < data[i])
    {
        max = data[i];
        index = idx[i];
    }
  }

  // all send my local max to master node
  // also send number of elements on this node, so master knows if I actually had data or not
  if(myRank == 0)
  {
    for(int i = 1; i < numPes; i++)
    {
        int64_t remotenp;
        int64_t remoteindex;
        double remotevalue;
        MPI_Recv(&remotenp, 1, MPI_INT64_T, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if(remotenp > 0)
        {
            // only do this if the remote process actually has data
            MPI_Recv(&remoteindex, 1, MPI_INT64_T, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&remotevalue, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(max < remotevalue)
            {
                // this is a new max value
                max = remotevalue;
                index = remoteindex;
            }
        }
    }
  }
  else
  {
      MPI_Send(&np, 1, MPI_INT64_T, 0, 0, MPI_COMM_WORLD);
      if(np > 0)
      {
        MPI_Send(&index, 1, MPI_INT64_T, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&max, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      }
  }

  // max and index will contain the overall max information
  // send to all ranks
  MPI_Bcast(&index, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return std::pair<int64_t, double>(index, max);
}


} // close namespace panache

