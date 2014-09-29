/*! \file
 * \brief Density fitting tensor generation and manipulation (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#include <fstream>
#include <algorithm>
#include <iostream>
#include "panache/ThreeIndexTensor.h"
#include "panache/FittingMetric.h"
#include "panache/Molecule.h"
#include "panache/BasisSet.h"
#include "panache/Lapack.h"
#include "panache/Exception.h"

// for reordering
#include "panache/MemorySwapper.h"
#include "panache/CNorm.h"
#include "panache/Orderings.h"

#include "panache/Output.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace panache
{

ThreeIndexTensor::ThreeIndexTensor(SharedBasisSet primary,
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


int ThreeIndexTensor::SetNThread(int nthread)
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



ThreeIndexTensor::~ThreeIndexTensor()
{
}

void ThreeIndexTensor::GenQso(int qflags, int storeflags)
{
    if(qso_)
        return; // already created!

    qso_ = StoredQTensorFactory(naux_, nso_, nso_, storeflags, "qso");
    qso_->Init();

    if(qflags & QGEN_DFQSO) 
        qso_->GenDFQso(fittingmetric_, primary_, auxiliary_, nthreads_);
//    else if(qflags & QGEN_CHQSO)
//        qso_->GenCHQso(fittingmetric_, primary_, auxiliary_, nthreads_);
    else
        throw RuntimeError("Error - unknown QSO type!");
}

void ThreeIndexTensor::SetCMatrix(double * cmo, int nmo, bool cmo_is_trans,
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
        reorder::CNorm * cnorm;

        if (order == BSORDER_GAMESS)
        {
            ord = new reorder::GAMESS_Ordering();
            cnorm = new reorder::GAMESS_CNorm();
        }
        else
            throw RuntimeError("Unknown ordering!");

        //std::cout << "BEFORE REORDERING:\n";
        //for(int i = 0; i < nmo_*nso_; i++)
        //    std::cout << Cmo_[i] << "\n";
        ReorderCMat(*ord, *cnorm);
        //std::cout << "AFTER REORDERING:\n";
        //for(int i = 0; i < nmo_*nso_; i++)
        //    std::cout << Cmo_[i] << "\n";
        delete ord;
    }
}


void ThreeIndexTensor::GenQTensors(int qflags, int storeflags)
{

#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
#endif

    // Make sure one of QGEN_CHQSO or QGEN_CHQSO is specified, but not both
    if(((qflags & QGEN_DFQSO) && (qflags & QGEN_CHQSO))
        || !((qflags & QGEN_DFQSO) || (qflags & QGEN_CHQSO)))
        throw RuntimeError("Set one and only one of QGEN_DFQSO or QGEN_CHQSO!");

    // Remove QSTORAGE_PACKED if specified
    storeflags &= ~(QSTORAGE_PACKED);

    // Always gen qso as packed and by q
    GenQso(qflags, storeflags | QSTORAGE_PACKED | QSTORAGE_BYQ);


    if((qflags & (15)) == 0)
        return; //?

    if(!Cmo_)
        throw RuntimeError("Set the c-matrix first!");

    std::vector<StoredQTensor::TransformMat> lefts;
    std::vector<StoredQTensor::TransformMat> rights;
    std::vector<StoredQTensor *> qouts;

    if(qflags & QGEN_QMO)
    {
        // generate Qmo
        qmo_ = StoredQTensorFactory(naux_, nmo_, nmo_, storeflags | QSTORAGE_PACKED, "qmo");
        qmo_->Init();
        qouts.push_back(qmo_.get());
        lefts.push_back(StoredQTensor::TransformMat(Cmo_.get(), nmo_));
        rights.push_back(StoredQTensor::TransformMat(Cmo_.get(), nmo_));
    }
    if(qflags & QGEN_QOO)
    {
        // generate Qoo
        qoo_ = StoredQTensorFactory(naux_, nocc_, nocc_, storeflags | QSTORAGE_PACKED, "qoo");
        qoo_->Init();
        qouts.push_back(qoo_.get());
        lefts.push_back(StoredQTensor::TransformMat(Cmo_occ_.get(), nocc_));
        rights.push_back(StoredQTensor::TransformMat(Cmo_occ_.get(), nocc_));
    }
    if(qflags & QGEN_QOV)
    {
        // generate Qov
        qov_ = StoredQTensorFactory(naux_, nocc_, nvir_, storeflags, "qov");
        qov_->Init();
        qouts.push_back(qov_.get());
        lefts.push_back(StoredQTensor::TransformMat(Cmo_occ_.get(), nocc_));
        rights.push_back(StoredQTensor::TransformMat(Cmo_vir_.get(), nvir_));
    }
    if(qflags & QGEN_QVV)
    {
        // generate Qvv
        qvv_ = StoredQTensorFactory(naux_, nvir_, nvir_, storeflags | QSTORAGE_PACKED, "qvv");
        qvv_->Init();
        qouts.push_back(qvv_.get());
        lefts.push_back(StoredQTensor::TransformMat(Cmo_vir_.get(), nvir_));
        rights.push_back(StoredQTensor::TransformMat(Cmo_vir_.get(), nvir_));
    }

    // actually do the transformation
    qso_->Transform(lefts, rights, qouts, nthreads_);

#ifdef PANACHE_TIMING
    tim.Stop();
    timer_genqtensors_.AddTime(tim);
#endif

}


int ThreeIndexTensor::GetQBatch_Base(double * outbuf, int bufsize, int qstart,
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


int ThreeIndexTensor::GetBatch_Base(double * outbuf, int bufsize, int ijstart,
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


std::unique_ptr<ThreeIndexTensor::StoredQTensor> & ThreeIndexTensor::ResolveTensorFlag(int tensorflag)
{
    switch(tensorflag)
    {
        case QGEN_DFQSO:
            return qso_;
        case QGEN_CHQSO:  // Stored in the same place as DFQSO
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


int ThreeIndexTensor::GetQBatch(int tensorflag, double * outbuf, int bufsize, int qstart)
{
    return GetQBatch_Base(outbuf, bufsize, qstart, ResolveTensorFlag(tensorflag).get());
}

int ThreeIndexTensor::GetBatch(int tensorflag, double * outbuf, int bufsize, int ijstart)
{
    return GetBatch_Base(outbuf, bufsize, ijstart, ResolveTensorFlag(tensorflag).get());
}


int ThreeIndexTensor::GetQBatch(int tensorflag, double * outbuf, int bufsize, QIterator qstart)
{
    if(qstart)
        return GetQBatch_Base(outbuf, bufsize, qstart.Index(), ResolveTensorFlag(tensorflag).get());
    else
        return 0;
}

int ThreeIndexTensor::GetBatch(int tensorflag, double * outbuf, int bufsize, IJIterator ijstart)
{
    if(ijstart)
        return GetBatch_Base(outbuf, bufsize, ijstart.Index(), ResolveTensorFlag(tensorflag).get());
    else
        return 0;
}


int ThreeIndexTensor::QBatchSize(int tensorflag)
{
    return ResolveTensorFlag(tensorflag)->ndim12();
}


int ThreeIndexTensor::BatchSize(int tensorflag)
{
    return naux_;
}

int ThreeIndexTensor::TensorDimensions(int tensorflag, int & naux, int & ndim1, int & ndim2)
{
    auto & qt = ResolveTensorFlag(tensorflag);
    naux = qt->naux();
    ndim1 = qt->ndim1();
    ndim2 = qt->ndim2();
    return qt->naux() * qt->ndim12();
}


bool ThreeIndexTensor::IsPacked(int tensorflag)
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

void ThreeIndexTensor::ReorderCMat(const reorder::Orderings & order, const reorder::CNorm & cnorm)
{
    using namespace reorder;
    using std::placeholders::_1;
    using std::placeholders::_2;

    // First, renormalize
    for(int i = 0; i < primary_->nshell(); i++)
    {
        const GaussianShell & s = primary_->shell(i);
        if(cnorm.NeedsCNorm(s.is_pure(), s.am()))
        {
            auto normfac = cnorm.GetCNorm(s.is_pure(), s.am());
            
            // multiply rows by factors
            for(size_t j = 0; j < normfac.size(); j++)
            {
                int jstart = s.function_index();

                for(int k = 0; k < nmo_; k++)
                    Cmo_[(jstart+j)*nmo_+k] *= normfac[j];
            }
        }
    }


    TotalMemorySwapper sf1(nmo_);  // swaps rows

    std::vector<PointerMap> vpm;

    //go through what would need to be changed in the primary basis
    for(int i = 0; i < primary_->nshell(); i++)
    {
        const GaussianShell & s = primary_->shell(i);
        if(order.NeedsInvReordering(s.is_pure(), s.am()))
            vpm.push_back(PointerMap(s.function_index(), order.GetInvOrder(s.is_pure(), s.am())));
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

void ThreeIndexTensor::SplitCMat(void)
{
    Cmo_occ_ = std::unique_ptr<double[]>(new double[nso_*nocc_]);
    Cmo_vir_ = std::unique_ptr<double[]>(new double[nso_*nvir_]);

    //std::fill(Cmo_occ_.get(), Cmo_occ_.get() + nso_*nocc_, 0.0);
    //std::fill(Cmo_vir_.get(), Cmo_vir_.get() + nso_*nvir_, 0.0);

    // note - Cmo_occ_ and Cmo_vir_ will always be in column major order!
    // Cmo_ is nso * nmo
    //! \todo BLAS call?
    for(int i = 0; i < nso_; i++)
    {
        // remember to remove the frozen orbitals
        for(int j = 0; j < nocc_; j++)
            Cmo_occ_[i*nocc_ + j] = Cmo_[i*nmo_+(j+nfroz_)];
        for(int j = 0; j < nvir_; j++)
            Cmo_vir_[i*nvir_ + j] = Cmo_[i*nmo_+(j+nocc_+nfroz_)];
    }
}

void ThreeIndexTensor::SetNOcc(int nocc, int nfroz)
{
    if(nocc <= 0)
        throw RuntimeError("Error - nocc <= 0!");

    if(Cmo_ == nullptr)
        throw RuntimeError("Error - C Matrix not set!");

    nocc_ = nocc;
    nfroz_ = nfroz;
    nvir_ = nmo_ - nocc - nfroz;

    SplitCMat();
}


void ThreeIndexTensor::PrintTimer(const char * name, const std::unique_ptr<StoredQTensor> & q) const
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

void ThreeIndexTensor::PrintTimings(void) const
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
    output::printf(std::string(93, '-').c_str());
    output::printf("\n\n");

    #endif
}


ThreeIndexTensor::IteratedQTensorByQ ThreeIndexTensor::IterateByQ(int tensorflag, double * buf, int bufsize)
{
    using std::placeholders::_1;

    // ugly
    IteratedQTensorByQ::GetBatchFunc gbf(std::bind(
                                   static_cast<int(ThreeIndexTensor::*)(int, double *, int, int)>(&ThreeIndexTensor::GetQBatch),
                                   this, tensorflag, buf, bufsize, _1));

    int ndim1, ndim2, naux;
    TensorDimensions(tensorflag, naux, ndim1, ndim2);

    IteratedQTensorByQ iqt(tensorflag, buf, bufsize, 
                           QBatchSize(tensorflag),
                           QIterator(naux, IsPacked(tensorflag)),
                           gbf);

    return iqt;
}


ThreeIndexTensor::IteratedQTensorByIJ ThreeIndexTensor::IterateByIJ(int tensorflag, double * buf, int bufsize)
{
    using std::placeholders::_1;

    // ugly
    IteratedQTensorByIJ::GetBatchFunc gbf(std::bind(
                                   static_cast<int(ThreeIndexTensor::*)(int, double *, int, int)>(&ThreeIndexTensor::GetBatch),
                                   this, tensorflag, buf, bufsize, _1));

    int ndim1, ndim2, naux;
    TensorDimensions(tensorflag, naux, ndim1, ndim2);

    IteratedQTensorByIJ iqt(tensorflag, buf, bufsize, 
                            BatchSize(tensorflag),
                            IJIterator(ndim1, ndim2, IsPacked(tensorflag)),
                            gbf);

    return iqt;
}


} // close namespace panache


