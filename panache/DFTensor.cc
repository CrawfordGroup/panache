/*! \file
 * \brief Density fitting tensor generation and manipulation (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#include <fstream>
#include <algorithm>
#include <iostream>
#include "panache/DFTensor.h"
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

DFTensor::DFTensor(SharedBasisSet primary,
                     SharedBasisSet auxiliary,
                     const std::string & directory,
                     int nthreads)
    : primary_(primary), auxiliary_(auxiliary), nthreads_(nthreads), directory_(directory)
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

    //remove trailing slashes
    while(directory_.size() > 1 && directory_.back() == '/')
        directory_ = directory_.substr(0, directory_.size()-1);

    Cmo_ = nullptr;
    nmo_ = 0;
    nmo2_ = 0;
    nso_ = primary_->nbf();
    nso2_ = nso_*nso_;
    naux_ = auxiliary_->nbf();
    nsotri_ = (nso_*(nso_+1))/2;
}


int DFTensor::SetNThread(int nthread)
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



DFTensor::~DFTensor()
{
}

void DFTensor::GenQso(int storeflags)
{
    if(qso_)
        return; // already created!

    qso_ = StoredQTensorFactory(naux_, nso_, nso_, storeflags, "qso");
    qso_->Init();

#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
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
            C_DGEMM('T','T',nm*nn, naux_, naux_, 1.0, B[threadnum], nm*nn, J, naux_, 0.0,
                    A[threadnum], naux_);


            // write to disk or store in memory
            for (int m0 = 0, m = mstart; m < mend; m0++, m++)
                    qso_->Write(A[threadnum] + (m0*nn)*naux_, nn, qso_->calcindex(m, nstart));
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
        C_DGEMM('T','T',nm*nm, naux_, naux_, 1.0, B[threadnum], nm*nm, J, naux_, 0.0,
                A[threadnum], naux_);


        // write to disk or store in memory
        /*
        for (int m0 = 0, m = mstart; m < mend; m0++, m++)
            for (int n0 = 0, n = mstart; n < mend; n0++, n++)
                if(m <= n)
                    qso_->Write(A[threadnum] + (m0*nm+n0)*naux_, 1, qso_->calcindex(m, n));
        */
        int iwrite = 1;
        for (int m0 = 0, m = mstart; m < mend; m0++, m++)
            qso_->Write(A[threadnum] + (m0*nm)*naux_, iwrite++, qso_->calcindex(m, mstart));

    }

    for(int i = 0; i < nthreads_; i++)
    {
        delete [] A[i];
        delete [] B[i];
    }

#ifdef PANACHE_TIMING
    tim.Stop();
    qso_->GenTimer().AddTime(tim);
#endif

}

void DFTensor::SetCMatrix(double * cmo, int nmo, bool cmo_is_trans,
                           int order)
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

    if(order != BSORDER_PSI4)
    {
        reorder::Orderings * ord;

        if (order == BSORDER_GAMESS)
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


void DFTensor::TransformQTensor(double * left, int lncols,
                                 double * right, int rncols,
                                 std::unique_ptr<StoredQTensor> & qout,
                                 int q,
                                 double * qso, double * qc, double * cqc)
{
#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
#endif
    // transform

    C_DGEMM('T', 'N', lncols, nso_, nso_, 1.0, left, lncols, qso, nso_, 0.0, qc, nso_);
    C_DGEMM('N', 'N', lncols, rncols, nso_, 1.0, qc, nso_, right, rncols, 0.0, cqc, rncols);

    // Write to memory or disk
    qout->WriteByQ(cqc, 1, q);

#ifdef PANACHE_TIMING
    tim.Stop();
    qout->GenTimer().AddTime(tim);
#endif
}


void DFTensor::GenQTensors(int qflags, int storeflags)
{
#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
#endif

    // temporary space
    // note that nmo, occupied, etc, are all <= nso
    // so this is the max space that would be needed
    double * qp = new double[nsotri_*nthreads_];   // packed q
    double * qe = new double[nso2_*nthreads_];   // expanded q
    double * qc = new double[nso2_*nthreads_];     // C(t) Q
    double * cqc = new double[nso2_*nthreads_];    // C(t) Q C

    // Remove QSTORAGE_PACKED if specified
    storeflags &= ~(QSTORAGE_PACKED);


    // Always gen qso as packed and by q
    GenQso(storeflags | QSTORAGE_PACKED | QSTORAGE_BYQ);


    if((qflags & (15)) == 0)
        return; //?

    if(!Cmo_)
        throw RuntimeError("Set the c-matrix first!");

    if(qflags & QGEN_QMO)
    {
        // generate Qmo
        qmo_ = StoredQTensorFactory(naux_, nmo_, nmo_, storeflags, "qmo");
        qmo_->Init();
    }
    if(qflags & QGEN_QOO)
    {
        // generate Qoo
        qoo_ = StoredQTensorFactory(naux_, nocc_, nocc_, storeflags, "qoo");
        qoo_->Init();
    }
    if(qflags & QGEN_QOV)
    {
        // generate Qov
        qov_ = StoredQTensorFactory(naux_, nocc_, nvir_, storeflags, "qov");
        qov_->Init();
    }
    if(qflags & QGEN_QVV)
    {
        // generate Qvv
        qvv_ = StoredQTensorFactory(naux_, nvir_, nvir_, storeflags, "qvv");
        qvv_->Init();
    }


    // Read in Qso, and create any requested tensors
#ifdef _OPENMP
    #pragma omp parallel for num_threads(nthreads_)
#endif
    for(int q = 0; q < naux_; q++)
    {
        int threadnum = 0;

#ifdef _OPENMP
        threadnum = omp_get_thread_num();
#endif

        double * myqp = qp + threadnum*nsotri_;
        double * myqe = qe + threadnum*nso2_;
        double * myqc = qc + threadnum*nso2_;
        double * mycqc = cqc + threadnum*nso2_;

        qso_->ReadByQ(myqp, 1, q);

        // expand packed matrix
        int index = 0;
        for(int i = 0; i < nso_; i++)
            for(int j = 0; j <= i; j++)
                myqe[i*nso_+j] = myqe[j*nso_+i] = myqp[index++];


        if(qflags & QGEN_QMO)
            TransformQTensor(Cmo_.get(), nmo_, Cmo_.get(), nmo_, qmo_, q, myqe, myqc, mycqc);
        if(qflags & QGEN_QOO)
            TransformQTensor(Cmo_occ_.get(), nocc_, Cmo_occ_.get(), nocc_, qoo_, q, myqe, myqc, mycqc);
        if(qflags & QGEN_QOV)
            TransformQTensor(Cmo_occ_.get(), nocc_, Cmo_vir_.get(), nvir_, qov_, q, myqe, myqc, mycqc);
        if(qflags & QGEN_QVV)
            TransformQTensor(Cmo_vir_.get(), nvir_, Cmo_vir_.get(), nvir_, qvv_, q, myqe, myqc, mycqc);
    }


    delete [] qp;
    delete [] qe;
    delete [] qc;
    delete [] cqc;

#ifdef PANACHE_TIMING
    tim.Stop();
    timer_genqtensors_.AddTime(tim);
#endif
}


int DFTensor::GetQBatch_Base(double * outbuf, int bufsize, int qstart,
                             StoredQTensor * qt)
{
#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
#endif

    int nq = (bufsize / qt->ndim12());

    if(nq == 0)
        throw RuntimeError("Error - buffer is to small to hold even one batch!");

    // get a batch
    int gotten = qt->ReadByQ(outbuf, nq, qstart);

#ifdef PANACHE_TIMING
    tim.Stop();
    qt->GetQBatchTimer().AddTime(tim);
#endif

    return gotten;
}


int DFTensor::GetBatch_Base(double * outbuf, int bufsize, int ijstart,
                             StoredQTensor * qt)
{
#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
#endif

    int nij = (bufsize / qt->naux());

    if(nij == 0)
        throw RuntimeError("Error - buffer is to small to hold even one batch!");

    // get a batch
    int gotten = qt->Read(outbuf, nij, ijstart);

#ifdef PANACHE_TIMING
    tim.Stop();
    qt->GetBatchTimer().AddTime(tim);
#endif

    return gotten;
}


std::unique_ptr<DFTensor::StoredQTensor> & DFTensor::ResolveTensorFlag(int tensorflag)
{
    switch(tensorflag)
    {
        case QGEN_QSO:
            return qso_;
        case QGEN_QMO:
            return qmo_;
        case QGEN_QOO:
            return qoo_;
        case QGEN_QOV:
            return qov_;
        case QGEN_QVV:
            return qvv_;
        default:
            throw RuntimeError("Unknown tensorflag");
    }
}


int DFTensor::GetQBatch(int tensorflag, double * outbuf, int bufsize, int qstart)
{
    return GetQBatch_Base(outbuf, bufsize, qstart, ResolveTensorFlag(tensorflag).get());
}

int DFTensor::GetBatch(int tensorflag, double * outbuf, int bufsize, int ijstart)
{
    return GetBatch_Base(outbuf, bufsize, ijstart, ResolveTensorFlag(tensorflag).get());
}


int DFTensor::QBatchSize(int tensorflag)
{
    return ResolveTensorFlag(tensorflag)->ndim12();
}


int DFTensor::BatchSize(int tensorflag)
{
    return naux_;
}

int DFTensor::TensorDimensions(int tensorflag, int & naux, int & ndim1, int & ndim2)
{
    auto & qt = ResolveTensorFlag(tensorflag);
    naux = qt->naux();
    ndim1 = qt->ndim1();
    ndim2 = qt->ndim2();
    return qt->naux() * qt->ndim12();
}


bool DFTensor::IsPacked(int tensorflag)
{
    return ResolveTensorFlag(tensorflag)->packed();
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

void DFTensor::ReorderCMat(reorder::Orderings & order)
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

void DFTensor::SplitCMat(void)
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

void DFTensor::SetNOcc(int nocc, int nfroz)
{
    if(nocc <= 0)
        throw RuntimeError("Error - nocc <= 0!");

    if(Cmo_ == nullptr)
        throw RuntimeError("Error - C Matrix not set!");

    nocc_ = nocc;
    nfroz_ = nfroz;
    nvir_ = nmo_ - nocc;

    SplitCMat();
}


void DFTensor::PrintTimer(const char * name, const std::unique_ptr<StoredQTensor> & q) const
{
    if(q)
    {
        // need to convert from std::atomic
        unsigned long genm = q->GenTimer().microseconds;
        unsigned long gent = q->GenTimer().timescalled;
        unsigned long getm = q->GetBatchTimer().microseconds;
        unsigned long gett = q->GetBatchTimer().timescalled;
        unsigned long getqm = q->GetQBatchTimer().microseconds;
        unsigned long getqt = q->GetQBatchTimer().timescalled;

        output::printf("%-6s  %17lu (%7lu)  %17lu (%7lu)  %17lu (%7lu)\n", name, 
                        genm, gent,
                        getm, gett,
                        getqm, getqt);
    }
    else
        output::printf("%-6s  %17s  %17s  %17s\n", name, "N/A", "N/A", "N/A"); 
}

void DFTensor::PrintTimings(void) const
{
    #ifdef PANACHE_TIMING

    output::printf("\n\n  ==> LibPANACHE DF Tensor Timings <==\n\n");
    output::printf("*Time (in microseconds), followed by number of times called in parentheses*\n");
    //output::printf(std::string(80, '-').c_str());
    output::printf("\n");
    output::printf("%-6s  %-27s  %-27s  %-27s\n",
                   "Tensor",
                   "        Generation", 
                   "         GetBatch", 
                   "         GetQBatch");
    output::printf("%-6s  %-27s  %-27s  %-27s\n", "------",
                   std::string(27,'-').c_str(), 
                   std::string(27,'-').c_str(), 
                   std::string(27,'-').c_str()); 
    PrintTimer("QSO", qso_);
    PrintTimer("QMO", qmo_);
    PrintTimer("QOO", qoo_);
    PrintTimer("QOV", qov_);
    PrintTimer("QVV", qvv_);
    output::printf(std::string(80, '-').c_str());
    output::printf("\n\n");

    #endif
}


} // close namespace panache


