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
#include "panache/Flags.h"
#include "panache/Iterator.h"
 
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


CyclopsQTensor::CyclopsQTensor(int storeflags, const std::string & name) : StoredQTensor(storeflags, name)
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

    // Note - storing this info in the class members, although
    // they are really only used in this function
    // (just to keep it somewhat consistent with GenCHQso)

    // Get the begin and end shells for this rank
    // Parallelization is over primary basis
    // but ending at shell boundaries
    auto shellrangeinfo = ShellRange2_(primary);

    mynelements_ = shellrangeinfo.first;
    myrange_ = shellrangeinfo.second;

    // The above calculation was only on the ij pair.
    // Now we include naux
    int64_t nelements = mynelements_ * naux;

    mydata_ = std::unique_ptr<double[]>(new double[nelements]);
    myidx_ = std::unique_ptr<int64_t[]>(new int64_t[nelements]);

    //std::cout << rank << ": MYRANGE: [" << range.first << " , " << range.second << ") NSHELL=" << primary->nshell() << "\n";
    //std::cout << rank << ": NELEMENTS: " << nelements << "\n";

    int64_t curidx = 0;

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
                    for (int m = 0; m < nm; m++)
                    for (int n = 0; n < nn; n++, index++)
                        B[threadnum][p*nm*nn + m*nn + n] = eribuffers[threadnum][index];
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
                    myidx_[curidx] = n*ndim1()*naux + m*naux + q;
                    mydata_[curidx] = A[threadnum][m0*naux*nn + n0*naux + q];
                    curidx++;
                }
            }
            else
            {
                for (int m0 = 0, m = mstart; m < mend; m0++, m++)
                for (int n0 = 0, n = nstart; n < nend; n0++, n++)
                for (int q = 0; q < naux; q++)
                {
                    myidx_[curidx] = n*ndim1()*naux + m*naux + q;
                    myidx_[curidx+1] = m*ndim1()*naux + n*naux + q;
                    mydata_[curidx] = mydata_[curidx+1] = A[threadnum][m0*naux*nn + n0*naux + q];
                    curidx += 2;
                }
            }
        }
    }

    //std::cout << rank << " CURIDX/NELEMENTS: " << curidx << " / " << nelements << "\n";

    tensor_->write(curidx, myidx_.get(), mydata_.get());

    for(int i = 0; i < nthreads; i++)
    {
        delete [] A[i];
        delete [] B[i];
    }

    mydata_.reset();
    myidx_.reset();

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
                    for (int on = 0; on <= om; on++)
                    {
                        myidx_[curidx] = ((om + mstart) * (om + mstart + 1))/2 + (on + nstart);
                        mydata_[curidx] = buffer[om * nN * nM * nN + on * nM * nN + om * nN + on];
                        curidx += 1;
                    }
                }
                else
                {    
                    for (int om = 0; om < nM; om++)
                    for (int on = 0; on < nN; on++)
                    {
                        myidx_[curidx] = ((om + mstart) * (om + mstart + 1))/2 + (on + nstart);
                        mydata_[curidx] = buffer[om * nN * nM * nN + on * nM * nN + om * nN + on];
                        curidx += 1;
                    }
                }
            }
        }
    }

    target.write(curidx, myidx_.get(), mydata_.get());
}

void CyclopsQTensor::ComputeRow_(std::vector<SharedTwoBodyAOInt> & eris, 
                                 int64_t row, CTF_Vector & target)
{
    SharedBasisSet basis = eris[0]->basis();

    auto ij = math::decomposeij_packed(row);
    int r = ij.first;
    int s = ij.second;
    int R = basis->function_to_shell(r);
    int S = basis->function_to_shell(s);

    int nR = basis->shell(R).nfunction();
    int nS = basis->shell(S).nfunction();
    int rstart = basis->shell(R).function_index();
    int sstart = basis->shell(S).function_index();

    int oR = r - rstart;
    int oS = s - sstart;

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
            int nint = integral->compute_shell(M,N,R,S);

            if(nint)
            {
                int nN = basis->shell(N).nfunction();
                int nstart = basis->shell(N).function_index();
    
                if (N == M)
                {
                    for (int om = 0; om < nM; om++)
                    for (int on = 0; on <= om; on++)
                    {
                        myidx_[curidx] = ((om + mstart) * (om + mstart + 1))/2 + (on + nstart);
                        mydata_[curidx] = buffer[om * nN * nR * nS + on * nR * nS + oR * nS + oS];
                        curidx += 1;
                    }
                }
                else
                {    
                    for (int om = 0; om < nM; om++)
                    for (int on = 0; on < nN; on++)
                    {
                        myidx_[curidx] = ((om + mstart) * (om + mstart + 1))/2 + (on + nstart);
                        mydata_[curidx] = buffer[om * nN * nR * nS + on * nR * nS + oR * nS + oS];
                        curidx += 1;
                    }
                }
            }
        }
    }

    target.write(curidx, myidx_.get(), mydata_.get());
}

void CyclopsQTensor::GenCHQso_(const SharedBasisSet primary,
                                                 double delta,
                                                 int nthreads)
{
    auto shellrangeinfo = ShellRange2_(primary);

    mynelements_ = shellrangeinfo.first;
    myrange_ = shellrangeinfo.second;

    // buffers are used in ComputeDiagonal_ and ComputeRow_
    mydata_ = std::unique_ptr<double[]>(new double[mynelements_]);
    myidx_ = std::unique_ptr<int64_t[]>(new int64_t[mynelements_]);

    std::vector<SharedTwoBodyAOInt> eris;

    for(int i = 0; i < nthreads; i++)
        eris.push_back(GetERI(primary, primary, primary, primary));

    int nQ = 0;
    int n = primary->nbf();
    int n12 = (n*(n+1))/2;


    CTF_Vector diag(n12, parallel::CTFWorld(), "CHODIAG");
    ComputeDiagonal_(eris, diag);

    // Temporary cholesky rows
    std::vector<CTF_Vector *> L;

    // List of selected pivots
    std::vector<int64_t> pivots;

    int64_t oneidx = 0;
    double oneval = 0;
    bool ismaster = parallel::IsMaster();

    while(nQ < n12)
    {
        auto maxel = FindVecMax_(diag);
        int64_t pivot = maxel.first;
        double Dmax = maxel.second;

        if(Dmax < delta || Dmax < 0.0) break;

        pivots.push_back(pivot); // should be ok distributed. They all have the same index/value

        double L_QQ = sqrt(Dmax);

        // no need to zero like in LocalQTensor - cyclops does it
        L.push_back(new CTF_Vector(n12, parallel::CTFWorld()));

        ComputeRow_(eris, pivot, *(L[nQ]));

        // [(m|Q) - L_m^P L_Q^P]
        oneidx = pivot;
        for(int P = 0; P < nQ; P++)
        {
            L[P]->read(1, &oneidx, &oneval);
            (*L[nQ])["i"] += (-1.0*oneval) * (*L[P])["i"];
        }

        // 1/L_QQ [(m|Q) - L_m^P L_Q^P]
        (*L[nQ]).scale(1.0 / L_QQ, "i");

        // Update the Schur complement diagonal
        diag["i"] -= (*L[nQ])["i"] * (*L[nQ])["i"];

        // for zeroing some elements
        std::vector<double> zeroval(nQ+1, 0.0);

        // pivot factor
        oneidx = pivot;
        oneval = L_QQ;

        if(ismaster)
        {
            // Zero the upper triangle
            L[nQ]->write(nQ+1, pivots.data(), zeroval.data());

            // Set the pivot factor
            L[nQ]->write(1, &oneidx, &oneval);

            // Force truly zero elements to zero
            diag.write(nQ+1, pivots.data(), zeroval.data());
        }
        else
        {
            // These are dummy writes. Only the master does writes
            // If not, value gets doubled (????)
            L[nQ]->write(0, pivots.data(), zeroval.data());
            L[nQ]->write(0, &oneidx, &oneval);
            diag.write(0, pivots.data(), zeroval.data());
        }

        nQ++;
    }

    // done with buffers
    mydata_.reset();
    myidx_.reset();
    eris.clear();
    pivots.clear();

    // Now we know the dimensions. initialize the tensor
    StoredQTensor::Init(nQ, n, n);

    // unfortunately we need a mapping of packed indices to unpacked
    // (for distributing the data at the end)
    // (well, we don't NEED it, but should be faster than decomposing the indices all the time)
    std::vector<std::pair<int64_t, int64_t>> packmap;
    packmap.reserve(n12);

    IJIterator ijit(n, n, true);
    for(int i = 0; i < n12; i++)
    {
        packmap.push_back(std::pair<int64_t, int64_t>(ijit.i(), ijit.j()));
        ++ijit;
    }
    

    // convert vector of CTF_Vector to a proper 3-index tensor
    // with unpacking
    double * localvals;
    int64_t * localidx;
    int64_t np;

    for(int i = 0; i < nQ; i++)
    {
        // need local data, plus a copy of the indices
        L[i]->read_local(&np, &localidx, &localvals);
        int64_t * idxcopy = new int64_t[np];
        std::copy(localidx, localidx+np, idxcopy);

        // row in 3-index tensor corresponds to i
        // get row/col from map
        // and remembering that this is a symmetric matrix, change
        // the index in localidx and reuse the data array
        for(int64_t j = 0; j < np; j++)
        {
            //auto pij = packmap[idxcopy[j]];
            auto pij = packmap.at(idxcopy[j]);
            localidx[j] = i + nQ * (pij.first + pij.second * n);
        }
        tensor_->write(np, localidx, localvals);

        // repeat for the other symmetric part
        for(int64_t j = 0; j < np; j++)
        {
            //auto pij = packmap[idxcopy[j]];
            auto pij = packmap.at(idxcopy[j]);
            localidx[j] = i + nQ * (pij.second + pij.first * n);
        }
        tensor_->write(np, localidx, localvals);

        free(localidx);
        free(localvals);
        delete [] idxcopy;
        delete L[i];
    }
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

  double max = 0;
  int64_t index = 0;

  if(np)
  {
      max = data[0];
      index = idx[0];
      for(int64_t i = 1; i < np; i++)
      {
          if(max < data[i])
          {
              max = data[i];
              index = idx[i];
          }
      }
  }

  free(idx);
  free(data);

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

void CyclopsQTensor::Finalize_(void)
{
}

void CyclopsQTensor::NoFinalize_(void)
{
}

} // close namespace panache

