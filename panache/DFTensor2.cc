/*! \file
 * \brief Density fitting tensor generation and manipulation (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#include <fstream>
#include <algorithm>
#include <iostream>
#include "panache/DFTensor2.h"
#include "panache/FittingMetric.h"
#include "panache/Molecule.h"
#include "panache/BasisSet.h"
#include "panache/Lapack.h"
#include "panache/Exception.h"

// for reordering
#include "panache/Orderings.h"
#include "panache/MemorySwapper.h"

#include "panache/ERI.h"

#include "panache/ERDERI.h"
#include "panache/Output.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace panache
{

DFTensor2::DFTensor2(SharedBasisSet primary,
                   SharedBasisSet auxiliary,
                   const std::string & filename,
                   int nthreads)
    : primary_(primary), auxiliary_(auxiliary), filename_(filename), nthreads_(nthreads) 
{
    output::printf("  ==> LibPANACHE DF Tensor <==\n\n");

    output::printf(" => Primary Basis Set <= \n\n");
    primary_->print_detail();

    output::printf(" => Auxiliary Basis Set <= \n\n");
    auxiliary_->print_detail();

    #ifdef _OPENMP
      if(nthreads_ <= 0)
          nthreads_ = omp_get_max_threads();
    #else
      nthreads_ = 1;
    #endif

    fittingmetric_ = std::shared_ptr<FittingMetric>(new FittingMetric(auxiliary_, nthreads_));
    fittingmetric_->form_eig_inverse();


    curq_ = 0;
    Cmo_ = nullptr;
    nmo_ = 0;
    nmo2_ = 0;
    nso_ = primary_->nbf();
    nso2_ = nso_*nso_;
    nsotri_ = ( (nso_ * (nso_+1) ) )/2;
    naux_ = auxiliary_->nbf();

    outbuffer_ = nullptr;
    outbuffersize_ = 0; 



}


int DFTensor2::SetNThread(int nthread)
{
    #ifdef _OPENMP
      if(nthread <= 0) 
          nthreads_ = omp_get_max_threads();
      else
         nthreads_ = nthread;
    #else
      nthreads_ = 1;
    #endif

    return nthreads_;
}



DFTensor2::~DFTensor2()
{
    CloseFile();
}

int DFTensor2::QsoDimensions(int & naux, int & nso2)
{
    nso2 = nso2_;
    naux = naux_;
    return nso2*naux;
}

void DFTensor2::SetCMatrix(double * cmo, int nmo, bool cmo_is_trans, 
                          BSOrder order)
{
    nmo_ = nmo;
    nmo2_ = nmo*nmo;

    Cmo_ = std::unique_ptr<double[]>(new double[nmo_*nso_]);

    if(cmo_is_trans)
    {
        for(int i = 0; i < nso_; i++)
        for(int j = 0; j < nmo_; j++)
            Cmo_[i*nmo_+j] = cmo[j*nso_+i];
    }
    else
        std::copy(cmo, cmo+(nmo_*nso_), Cmo_.get());

    if(order != BSOrder::Psi4)
    {
        reorder::Orderings * ord;

        if (order == BSOrder::GAMESS)
            ord = new reorder::GAMESS_Ordering();
        else
            throw RuntimeError("Unknown ordering!");

        //std::cout << "BEFORE REORDERING:\n";
        //for(int i = 0; i < nmo_*nso_; i++)
        //    std::cout << Cmo_[i] << "\n";
        ReorderCMat(*ord);
        //std::cout << "AFTER REORDERING:\n";
        //for(int i = 0; i < nmo_*nso_; i++)
        //    std::cout << Cmo_[i] << "\n";
        delete ord;
    }    
}


void DFTensor2::GenQso(bool inmem)
{

#ifdef PANACHE_TIMING
    timer_genqso.Reset();
    timer_genqso.Start();
#endif

    int maxpershell = primary_->max_function_per_shell();
    int maxpershell2 = maxpershell*maxpershell;

    double * J = fittingmetric_->get_metric();

    // default constructor = zero basis
    SharedBasisSet zero(new BasisSet);

    std::vector<std::shared_ptr<TwoBodyAOInt>> eris;
    std::vector<const double *> eribuffers;
    std::vector<double *> A, B;

    for(int i = 0; i < nthreads_; i++)
    {
        eris.push_back(GetERI(auxiliary_,zero,primary_,primary_));
        eribuffers.push_back(eris[eris.size()-1]->buffer());

        // temporary buffers
        A.push_back(new double[naux_*maxpershell2]);
        B.push_back(new double[naux_*maxpershell2]);
    }

    isinmem_ = inmem; // store this so we know later

    qso_.reset();
    curq_ = 0;
    CloseFile(); // safe even if not opened

    // Allocate memory or open the file on disk
    if(!isinmem_)
    {
        OpenFile();
        ResetFile();
    }
    else
        qso_ = std::unique_ptr<double[]>(new double[nsotri_ * naux_]);




    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads_)
    #endif
    for (int M = 0; M < primary_->nshell(); M++)
    {
        int threadnum = 0;
        #ifdef _OPENMP
            threadnum = omp_get_thread_num();
        #endif

        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        int mend = mstart + nm;

        for (int N = 0; N < M; N++)
        {
            int nn = primary_->shell(N).nfunction();
            int nstart = primary_->shell(N).function_index();
            //int nend = nstart + nn;

            for (int P = 0; P < auxiliary_->nshell(); P++)
            {
                int np = auxiliary_->shell(P).nfunction();
                int pstart = auxiliary_->shell(P).function_index();
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
            C_DGEMM('N','N',naux_, nm*nn, naux_, 1.0, J, naux_, B[threadnum], nm*nn, 0.0,
                    A[threadnum], nm*nn);


            // write to disk or store in memory
            //! \todo rearrange to that writes are more sequential?
            if(isinmem_)
            {
                // safe to write to qso_ in parallel - each thread writes to a different spot
                for (int p = 0; p < naux_; p++)
                {
                    for (int m0 = 0, m = mstart; m < mend; m0++, m++)
                        std::copy(A[threadnum]+p*nm*nn+m0*nn,
                                  A[threadnum]+p*nm*nn+m0*nn+nn, 
                                  qso_.get() + p*nsotri_ + ((m*(m+1))>>1) + nstart);
                }
            }
            else
            {
                for (int p = 0; p < naux_; p++)
                {
                    for (int m0 = 0, m = mstart; m < mend; m0++, m++)
                    {
                        // must be done in a critical section
                        #ifdef _OPENMP
                        #pragma omp critical
                        #endif
                        {
                            matfile_->seekp(sizeof(double)*(p*nsotri_ + ((m*(m+1))>>1) + nstart), std::ios_base::beg);
                            matfile_->write(reinterpret_cast<const char *>(A[threadnum] + p*nm*nn + m0*nn), nn*sizeof(double));
                        }
                    }
                }
            }
        }


        // Special Case: N = M
        for (int P = 0; P < auxiliary_->nshell(); P++)
        {
            int np = auxiliary_->shell(P).nfunction();
            int pstart = auxiliary_->shell(P).function_index();
            int pend = pstart + np;

            eris[threadnum]->compute_shell(P,0,M,M);

            for (int p = pstart, index = 0; p < pend; p++)
            {
                for (int m = 0; m < nm; m++)
                {
                    for (int n = 0; n < nm; n++, index++)
                    {
                        B[threadnum][p*nm*nm + m*nm + n] = eribuffers[threadnum][index];
                    }
                }
            }
        }

        // we now have a set of columns of B, although "condensed"
        // we can do a DGEMM with J
        C_DGEMM('N','N',naux_, nm*nm, naux_, 1.0, J, naux_, B[threadnum], nm*nm, 0.0,
                A[threadnum], nm*nm);


        // write to disk
        //! \todo rearrange to that writes are more sequential?
        if(isinmem_)
        {
            // safe to write to qso_ in parallel - each thread writes to a different spot
            for (int p = 0; p < naux_; p++)
            {
                int nwrite = 1;
                for (int m0 = 0, m = mstart; m < mend; m0++, m++)
                {
                    std::copy(A[threadnum]+p*nm*nm+m0*nm,
                              A[threadnum]+p*nm*nm+m0*nm+nwrite,
                              qso_.get()+p*nsotri_ + ((m*(m+1))>>1) + mstart);
                    nwrite += 1;
                }
            }
        }
        else
        {
            for (int p = 0; p < naux_; p++)
            {
                int nwrite = 1;
                for (int m0 = 0, m = mstart; m < mend; m0++, m++)
                {
                    // must be done in a critical section
                    #ifdef _OPENMP
                    #pragma omp critical
                    #endif
                    {
                        matfile_->seekp(sizeof(double)*(p*nsotri_ + ((m*(m+1))>>1) + mstart), std::ios_base::beg);
                        matfile_->write(reinterpret_cast<const char *>(A[threadnum] + p*nm*nm + m0*nm), nwrite*sizeof(double));
                        nwrite += 1;
                    }
                }
            }
        }

    }

    if(!isinmem_)
        ResetFile();

    for(int i = 0; i < nthreads_; i++)
    {
        delete [] A[i];
        delete [] B[i];
    }

#ifdef PANACHE_TIMING
    timer_genqso.Stop();
    output::printf("  **TIMER: DFTensor2 Total GenQso (%s): %lu (%lu calls)\n",
                   (inmem ? "CORE" : "DISK"),
                   timer_genqso.Microseconds(),
                   timer_genqso.TimesCalled());

    // reset the other timers
    timer_getbatch_qso.Reset();
    timer_getbatch_qmo.Reset();
    timer_getbatch_qov.Reset();
    timer_getbatch_qoo.Reset();
    timer_getbatch_qvv.Reset();
#endif

}



void DFTensor2::SetOutputBuffer(double * buf, long int size)
{
    outbuffer_ = buf;
    outbuffersize_ = size;
}


#ifdef PANACHE_DISKPREFETCH
int DFTensor2::GetBatch_Base(int ntoget)
{
    int actualtoget = std::min(ntoget, naux_ - curq_);

    if(actualtoget <= 0)
        return 0;

    if(isinmem_)
    {
        std::copy(qso_.get() + curq_*nsotri_, qso_.get() + (curq_+actualtoget)*nsotri_, q_.get());
        curq_ += actualtoget;
    }
    else
    {
        if(curq_ == 0)
        {
            // first time through - just grab the data
            matfile_->read(reinterpret_cast<char *>(q_.get()), actualtoget*nsotri_*sizeof(double));
        }
        else
        {
            // wait for the filling of q2_ to finish
            int futureget = fill_future_.get();

            if(futureget != actualtoget)
            {
                std::cout << "Error: future returned " << futureget << " but I expected " << actualtoget << "\n";
                throw RuntimeError("Async error");
            }

            // swap
            std::swap(q_, q2_);
        }

        curq_ += actualtoget;

        // async get the next batch into q2
        int futuretoget = std::min(ntoget, naux_ - curq_);

        if(futuretoget != 0)
        {
            fill_future_ = std::async(std::launch::async,                                 // policy
                                      [] (int toget, size_t bytes, double * outbuf, std::fstream * file)  // function/lambda
                                         { 
                                            //std::cout << "PREFETCHING " << toget << " MATRICES ( " << bytes << " BYTES)...";
                                            file->read(reinterpret_cast<char *>(outbuf), bytes); 
                                            //std::cout << "Done\n";
                                            return toget;
                                         },
                                      futuretoget, futuretoget*nsotri_*sizeof(double), q2_.get(), matfile_.get()  // args to lambda
                                    );
        }
    }


    return actualtoget;
}

#else

int DFTensor2::GetBatch_Base(int ntoget)
{
    // all reads should be sequential, therefore we can just start at the beginning
    //  (reset file at end of GenQ) and go from there

    int actualtoget = std::min(ntoget, naux_ - curq_);

    if(actualtoget <= 0)
        return 0;

    if(isinmem_)
        std::copy(qso_.get() + curq_*nsotri_, qso_.get() + (curq_+actualtoget)*nsotri_, q_.get());
    else
        matfile_->read(reinterpret_cast<char *>(q_.get()), actualtoget*nsotri_*sizeof(double));

    curq_ += actualtoget;

    return actualtoget;
}
#endif



int DFTensor2::GetBatch_Qso(void)
{

#ifdef PANACHE_TIMING
    timer_getbatch_qso.Start();
#endif

    int nq = (outbuffersize_ / nso2_);

    // is this the first batch?
    if(curq_ == 0)
    {
        if(outbuffersize_ < nso2_)
            throw RuntimeError("Error - buffer is to small to hold even one row!");

        // allocate buffers for some q
        q_ = std::unique_ptr<double[]>(new double[nq * nsotri_]);

        #ifdef PANACHE_DISKPREFETCH
        q2_ = std::unique_ptr<double[]>(new double[nq * nsotri_]);
        #endif
    }


    // get a batch    
    int gotten;
    if((gotten = GetBatch_Base(nq)))
    {
        // expand into mat
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads_)
        #endif
        for(int i = 0; i < gotten; i++)
        {
            double * my_q = q_.get() + i*nsotri_;

            int index = 0;
            for(int j = 0; j < nso_; j++)
                for(int k = 0; k <= j; k++)
                    outbuffer_[i*nso2_ + k*nso_ + j] 
                             = outbuffer_[i*nso2_ + j*nso_ + k]
                             = my_q[index++];
        }
    }


#ifdef PANACHE_TIMING
    timer_getbatch_qso.Stop();
    if(gotten == 0)
    {
        output::printf("  **TIMER: DFTensor2 Total GetBatch_Qso (%s): %lu (%lu calls)\n",
                       (isinmem_ ? "CORE" : "DISK"),
                       timer_getbatch_qso.Microseconds(),
                       timer_getbatch_qso.TimesCalled());
    }
#endif

    if(gotten == 0)
    {
        // free memory
        q_.reset();

        #ifdef PANACHE_DISKPREFETCH
        q2_.reset();
        #endif
    }

    return gotten;
}



int DFTensor2::GetBatch_transform(double * left, int lncols, 
                                 double * right, int rncols,
                                 Timer & timer, const char * timername,
                                 int nthreads)
{
#ifdef PANACHE_TIMING
    timer.Start();
#endif

    int batchsize = lncols * rncols;
    int nq = (outbuffersize_ / batchsize );


    // first batch?
    if(curq_ == 0)
    {
        if(outbuffersize_ < batchsize)
            throw RuntimeError("Error - buffer is to small to hold even one row!");

        // allocate buffers for some q
        q_ = std::unique_ptr<double[]>(new double[nq * nsotri_]);

        #ifdef PANACHE_DISKPREFETCH
        q2_ = std::unique_ptr<double[]>(new double[nq * nsotri_]);
        #endif

        q_single_ = std::unique_ptr<double[]>(new double[nq * nso2_]);

        //  QC, then (Ct) QC, with the right part stored in qc_
        qc_ = std::unique_ptr<double[]>(new double[nq * rncols * nso_]);
    }


    int gotten;
    if((gotten = GetBatch_Base(nq)))
    {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for(int i = 0; i < gotten; i++)
        {
            double * my_q = q_.get() + i*nsotri_;
            double * my_q_single = q_single_.get() + i*nso2_;

            double * my_qc;

            my_qc = qc_.get() + i*nso_*rncols;

            // Apply the matrices
            // Keep in mind that for the first step q_ is a (lower) part of a symmetric matrix
            // so first expand into a buffer. But we only need to fill half and we can use DSYMM.
            // We will fill the lower triangle.
            int index = 0;
            for(int j = 0; j < nso_; j++)
                for(int k = 0; k <= j; k++)
                    my_q_single[j*nso_ + k] = my_q[index++];


            // matrix multiply w/ symmetric matrix
            C_DSYMM('L','L',nso_,rncols,1.0, my_q_single, nso_, right, rncols, 0.0, my_qc, rncols);

            // Then regular matrix multiplication
            C_DGEMM('T','N',lncols, rncols, nso_, 1.0, left, lncols, my_qc, rncols, 0.0, outbuffer_ + i*lncols*rncols, rncols);
        }
    }

    if(gotten == 0)
    {
        // free memory
        q_.reset();
        q_single_.reset();
        qc_.reset();

        #ifdef PANACHE_DISKPREFETCH
        q2_.reset();
        #endif

    }

#ifdef PANACHE_TIMING
    timer.Stop();
    if(gotten == 0)
    {
        output::printf("  **TIMER: DFTensor2 Total %s (%s): %lu (%lu calls)\n",
                       timername,
                       (isinmem_ ? "CORE" : "DISK"),
                       timer.Microseconds(),
                       timer.TimesCalled());
    }
#endif

    return gotten;
}




int DFTensor2::GetBatch_Qmo(void)
{

    if(!Cmo_)
        throw RuntimeError("Error - I don't have a C matrix!");

    int gotten = GetBatch_transform(Cmo_.get(), nmo_,
                                    Cmo_.get(), nmo_,
                                    timer_getbatch_qmo, "GetBatch_Qmo",
                                    nthreads_);

    return gotten;
}



int DFTensor2::GetBatch_Qov(void)
{
    if(!Cmo_)
        throw RuntimeError("Error - I don't have a C matrix!");
    if(!Cmo_occ_ || !Cmo_vir_)
        throw RuntimeError("Error - Set occupied and virtual orbitals first!");

    int gotten = GetBatch_transform(Cmo_occ_.get(), nocc_,
                                    Cmo_vir_.get(), nvir_,
                                    timer_getbatch_qov, "GetBatch_Qov",
                                    nthreads_);

    return gotten;
}

int DFTensor2::GetBatch_Qoo(void)
{
    if(!Cmo_)
        throw RuntimeError("Error - I don't have a C matrix!");
    if(!Cmo_occ_ || !Cmo_vir_)
        throw RuntimeError("Error - Set occupied and virtual orbitals first!");

    int gotten = GetBatch_transform(Cmo_occ_.get(), nocc_,
                                    Cmo_occ_.get(), nocc_,
                                    timer_getbatch_qoo, "GetBatch_Qoo",
                                    nthreads_);

    return gotten;
}

int DFTensor2::GetBatch_Qvv(void)
{
    if(!Cmo_)
        throw RuntimeError("Error - I don't have a C matrix!");
    if(!Cmo_occ_ || !Cmo_vir_)
        throw RuntimeError("Error - Set occupied and virtual orbitals first!");

    int gotten = GetBatch_transform(Cmo_vir_.get(), nvir_,
                                    Cmo_vir_.get(), nvir_,
                                    timer_getbatch_qvv, "GetBatch_Qvv",
                                    nthreads_);

    return gotten;
}


// note - passing by value for the vector
static void Reorder(std::vector<unsigned short> order, std::vector<double *> pointers,
                    reorder::MemorySwapper & sf)
{
    long int size = order.size();

    // original order is 1 2 3 4 5 6....
    std::vector<unsigned short> currentorder(size);

    for(long int i = 0; i < size; i++)
        currentorder[i] = i+1;

    for(long int i = 0; i < size; i++)
    {
        // find the index in the current order
        long int cindex = 0;
        bool found = false;

        for(int j = 0; j < size; j++)
        {
            if(currentorder[j] == order[i])
            {
                found = true;
                cindex = j;
                break;
            }
        }
        if(!found)
            throw RuntimeError("Error in reordering - index not found?");


        // we shouldn't swap anything that was previously put in place...
        if(cindex < i)
            throw RuntimeError("Error in reordering - going to swap something I shouldn't");

        //swap
        if(cindex != i)
        {
            sf.swap(pointers[i], pointers[cindex]);
            std::swap(currentorder[i], currentorder[cindex]);
        }
    }

    // double check
    for(long int i = 0; i < size; i++)
    {
        if(currentorder[i] != order[i])
            throw RuntimeError("Reordering failed!");
    }
}

void DFTensor2::ReorderCMat(reorder::Orderings & order)
{
    using namespace reorder;
    using std::placeholders::_1;
    using std::placeholders::_2;

    TotalMemorySwapper sf1(nmo_);  // swaps rows

    std::vector<PointerMap> vpm;

    //go through what would need to be changed in the primary basis
    for(int i = 0; i < primary_->nshell(); i++)
    {
        const GaussianShell & s = primary_->shell(i);
        if(s.is_pure())
        {
            if(order.NeedsInvSphReordering(s.am()))
                vpm.push_back(PointerMap(s.function_index(), order.GetInvSphOrder(s.am())));
        }
        else
        {
            if(order.NeedsInvCartReordering(s.am()))
                vpm.push_back(PointerMap(s.function_index(), order.GetInvCartOrder(s.am())));
        }
    }


    std::vector<double *> pointers(primary_->max_function_per_shell());


    // Swap rows
    for(auto & it : vpm)
    {
        size_t ntoswap = it.order.size();

        for(size_t n = 0; n < ntoswap; n++)
            pointers[n] = &(Cmo_[(it.start+n)*nmo_]);

        Reorder(it.order, pointers, sf1);
    }
}

void DFTensor2::SplitCMat(void)
{
    Cmo_occ_ = std::unique_ptr<double[]>(new double[nso_*nocc_]);
    Cmo_vir_ = std::unique_ptr<double[]>(new double[nso_*nvir_]);

    std::fill(Cmo_occ_.get(), Cmo_occ_.get() + nso_*nocc_, 0.0);
    std::fill(Cmo_vir_.get(), Cmo_vir_.get() + nso_*nvir_, 0.0);

    // note - Cmo_occ_ and Cmo_vir_ will always be in column major order!
    // Cmo_ is nso * nmo
    //! \todo BLAS call?
    for(int i = 0; i < nso_; i++)
    {
        for(int j = 0; j < nocc_; j++)
            Cmo_occ_[i*nocc_ + j] = Cmo_[i*nmo_+j];
        for(int j = 0; j < nvir_; j++)
            Cmo_vir_[i*nvir_ + j] = Cmo_[i*nmo_+(j+nocc_)];
    }
}

void DFTensor2::SetNOcc(int nocc)
{
    if(nocc <= 0)
        throw RuntimeError("Error - nocc <= 0!");

    if(Cmo_ == nullptr)
        throw RuntimeError("Error - C Matrix not set!");

    nocc_ = nocc;
    nvir_ = nmo_ - nocc;

    SplitCMat();
}

} // close namespace panache
