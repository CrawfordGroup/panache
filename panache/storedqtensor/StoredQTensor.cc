/*! \file
 * \brief Generic three-index tensor storage (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/storedqtensor/StoredQTensor.h"
#include "panache/Exception.h"
#include "panache/Flags.h"

namespace panache
{

StoredQTensor::StoredQTensor(int storeflags, const std::string & name)
               : name_(name), storeflags_(storeflags)
{
    filled_ = false;
}

int StoredQTensor::naux(void) const
{
    return naux_;
}

int StoredQTensor::ndim1(void) const
{
    return ndim1_;
}

int StoredQTensor::ndim2(void) const
{
    return ndim2_;
}

int StoredQTensor::ndim12(void) const
{
    return ndim12_;
}

int StoredQTensor::storesize(void) const
{
    return ndim12_*naux_;
}

int StoredQTensor::storeflags(void) const
{
    return storeflags_;
}

int StoredQTensor::packed(void) const
{
    return (storeflags_ & QSTORAGE_PACKED);
}

int StoredQTensor::byq(void) const
{
    return (storeflags_ & QSTORAGE_BYQ);
}

bool StoredQTensor::filled(void) const
{
    return filled_;
}

void StoredQTensor::markfilled(void)
{
    filled_ = true;
}

const std::string & StoredQTensor::name(void) const
{
    return name_;
}

int StoredQTensor::calcindex(int i, int j) const
{
    if(!packed())
        return (i*ndim2_+j);
    else if(i >= j)
        return ((i*(i+1))>>1) + j;
    else
        return ((j*(j+1))>>1) + i;
}

void StoredQTensor::Init(int naux, int ndim1, int ndim2)
{
    naux_ = naux;
    ndim1_ = ndim1;
    ndim2_ = ndim2;

    if(packed() && ndim1 != ndim2)
        throw RuntimeError("non square packed matrices?");

    ndim12_ = (packed() ? (ndim1_ * (ndim2_+1))/2 : ndim1_*ndim2_);

    Init_();
}


void StoredQTensor::Init(const StoredQTensor & rhs)
{
    Init(rhs.naux_, rhs.ndim1_, rhs.ndim2_);
}


StoredQTensor::~StoredQTensor()
{
}

int StoredQTensor::StoreFlags(void) const
{
    return storeflags_;
}

int StoredQTensor::Read(double * data, int nij, int ijstart)
{
    // note - don't check (nij > 0)!
    // it may occassionally be called with nij = 0

    if(nij < 0)
        throw RuntimeError("Read() passed with negative nij!");

    if(ijstart + nij >= ndim12())
        nij = ndim12() - ijstart;

    Read_(data, nij, ijstart);
    return nij;
}

int StoredQTensor::ReadByQ(double * data, int nq, int qstart)
{
    // note - don't check (nq > 0)!
    // it may occassionally be called with nq = 0

    if(nq < 0)
        throw RuntimeError("Read() passed with negative nq!");

    if(qstart + nq >= naux_)
        nq = naux_-qstart;

    ReadByQ_(data, nq, qstart);
    return nq;
}


CumulativeTime & StoredQTensor::GenTimer(void)
{
    return gen_timer_;
}

CumulativeTime & StoredQTensor::GetQBatchTimer(void)
{
    return getqbatch_timer_;
}

CumulativeTime & StoredQTensor::GetBatchTimer(void)
{
    return getijbatch_timer_;
}

void StoredQTensor::GenDFQso(const SharedFittingMetric fit,
                                     const SharedBasisSet primary,
                                     const SharedBasisSet auxiliary,
                                     int nthreads)
{
#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
#endif

    GenDFQso_(fit, primary, auxiliary, nthreads);
    filled_ = true;

#ifdef PANACHE_TIMING
    tim.Stop();
    GenTimer().AddTime(tim);
#endif
}

void StoredQTensor::GenCHQso(const SharedBasisSet primary,
                                     double delta,
                                     int nthreads)
{
#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
#endif

    GenCHQso_(primary, delta, nthreads);
    filled_ = true;

#ifdef PANACHE_TIMING
    tim.Stop();
    GenTimer().AddTime(tim);
#endif
}


void StoredQTensor::Transform(const std::vector<TransformMat> & left,
                                        const std::vector<TransformMat> & right,
                                        std::vector<StoredQTensor *> results,
                                        int nthreads)
{
    Transform_(left, right, results, nthreads);
    for(auto it : results)
        it->filled_ = true;
}

void StoredQTensor::Finalize(int nthreads)
{
// Finalizing is part of generation
#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
#endif

    // call derived class function
    Finalize_(nthreads);

#ifdef PANACHE_TIMING
    tim.Stop();
    GenTimer().AddTime(tim);
#endif
}

void StoredQTensor::NoFinalize(void)
{
    // call derived class function
    NoFinalize_();
}

} // close namespace panache

