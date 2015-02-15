/*! \file
 * \brief Three-index tensor storage in memory (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/storedqtensor/MemoryQTensor.h"
#include "panache/storedqtensor/DiskQTensor.h"
#include "panache/Lapack.h"
#include "panache/ERI.h"
#include "panache/Flags.h"
#include "panache/Exception.h"
#include "panache/BasisSet.h"
#include "panache/FittingMetric.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace panache
{


void MemoryQTensor::Write_(double * data, int nij, int ijstart)
{
    int inaux = naux();

    if(byq())
    {
        int indim12 = ndim12();

        for(int q = 0, qoff = 0; q < inaux; q++, qoff += indim12)
            for(int ij = 0, ijoff = 0; ij < nij; ij++, ijoff += inaux)
                data_[qoff+ijstart+ij] = data[ijoff+q];
    }
    else
    {
        double * start = data_.get() + ijstart * inaux;
        std::copy(data, data+inaux, start);
    }
}


void MemoryQTensor::WriteByQ_(double * data, int nq, int qstart, bool ijpacked)
{
    ////////////////////////////
    // wow is this thing a mess
    // hope you know index math
    ///////////////////////////

    //! \todo Know exactly where inaux, etc, are used
    const int indim12 = ndim12();
    int inaux = naux();

    if(ijpacked == packed())
    {
        if(byq())
        {
            std::copy(data,
                      data + nq*indim12,
                      data_.get()+qstart*indim12);
        }
        else
        {
            const int inaux = naux();
    
            for(int q = 0, qoff = 0; q < nq; q++, qoff += indim12)
                for(int ij = 0, ijoff = 0; ij < indim12; ij++, ijoff += inaux)
                    data_[ijoff+qstart+q] = data[qoff+ij];
        }
    }
    else if(ijpacked)
    {
        const int indim1 = ndim1();
        const int indim2 = ndim2();

        // ijpacked, but not packed()
        // unpack ij
        if(byq())
        {
            for(int q0 = 0, q = qstart; q0 < nq; q0++)
            {
                const int qq = q*indim12;
                const int qq0 = q0*indim12;

                for(int i = 0; i < indim1; i++)
                {
                    const int ii = i*indim2;

                    for(int j = 0; j <= i; j++)
                      data_[qq + ii + j] = data_[qq + j*indim1 + i] 
                                         = data[qq0 + (i*(i+1)<<1) + j];
                }
            }
        }
        else // not stored by q
        {
            const int indim12_tmp = (indim1*(indim1+1))/2;

            for(int i = 0; i < indim1; i++)
            {
                const int ii = i*indim2*inaux;
                const int iia = i*inaux;

                for(int j = 0; j <= i; j++)
                {
                    const int jj = j*indim1*inaux;
                    const int jja = j*inaux;

                    for(int q0 = 0, q = qstart; q0 < nq; q0++)
                        data_[ii + jja + q] = data_[jj + iia + q] 
                                            = data[q0*indim12_tmp + (i*(i+1)<<1) + j];
                }
            }
        }
    }
    else // not ijpacked, but packed()
    {
        const int indim1 = ndim1();
        const int indim2 = ndim2();

        // ijpacked, but not packed()
        // unpack ij
        if(byq())
        {
            for(int q0 = 0, q = qstart; q0 < nq; q0++)
            {
                const int qq = q*indim12;
                const int qq0 = q0*indim12;

                for(int i = 0; i < indim1; i++)
                {
                    //const int ii = i*indim2;
                    const int i12 = (i*(i+1))/2;

                    //! \todo could be replaced with std::copy
                    for(int j = 0; j <= i; j++)
                      data_[qq + i12 + j] = data[qq0 + i*indim2 + j];
                }
            }
        }
        else // not stored by q
        {
            for(int i = 0; i < indim1; i++)
            {
                const int i12 = (i*(i+2))/2;

                for(int j = 0; j <= i; j++)
                {
                    const int ij12 = (i12 + j)*inaux;

                    for(int q0 = 0, q = qstart; q0 < nq; q0++)
                        data_[ij12 + q] = data[q0*indim1*indim2 + i*indim2 + j];
                }
            }
        }
    }
}

void MemoryQTensor::Read_(double * data, int nij, int ijstart)
{
    int inaux = naux();

    if(byq())
    {
        int indim12 = ndim12();

        for(int q = 0, qoff = 0; q < inaux; q++, qoff += indim12)
            for(int ij = 0, ijoff = 0; ij < nij; ij++, ijoff += inaux)
                data[ijoff+q] = data_[qoff+ijstart+ij];
    }
    else
    {
        double * start = data_.get() + ijstart * inaux;
        std::copy(start, start+nij*inaux, data);
    }
}

void MemoryQTensor::ReadByQ_(double * data, int nq, int qstart)
{
    int indim12 = ndim12();

    if(byq())
    {
        double * start = data_.get()+qstart*indim12;
        std::copy(start, start+nq*indim12, data);
    }
    else
    {
        int inaux = naux();

        for(int q = 0, qoff = 0; q < nq; q++, qoff += indim12)
            for(int ij = 0, ijoff = 0; ij < indim12; ij++, ijoff += inaux)
                data[qoff+ij] = data_[ijoff+qstart+q];
    }
}


void MemoryQTensor::Init_(void)
{
    if(!data_)
        data_ = std::unique_ptr<double []>(new double[storesize()]);
}


MemoryQTensor::MemoryQTensor(int storeflags, const std::string & name, const std::string & directory) 
    : LocalQTensor(storeflags, name, directory)
{
    // does the file exist?
    // (set by LocalQTensor constructor)
    if(existed_ && (storeflags & QSTORAGE_READDISK))
    {
        // note - don't do diskqt(this) since this isn't completely constructed yet!
        // note2 - it's ok if the INMEM vs. DISK flags are incorrect in diskqt. They
        // aren't saved to the .dim file or matter much after construction
        DiskQTensor diskqt(storeflags, name, directory);

        if(!diskqt.filled())
            throw RuntimeError("Can't read from disk, but file exists?");

        // Init sizes, etc, 
        // from StoredQTensor base class
        // (calls MemoryQTensor::Init_)
        Init(diskqt);

        if(byq())
          diskqt.ReadByQ(data_.get(), naux(), 0);
        else
          diskqt.Read(data_.get(), ndim12(), 0);

        markfilled();

        // may delete file, based on store flags
    }
}

MemoryQTensor::~MemoryQTensor()
{
    if(!existed_ && (storeflags() & QSTORAGE_KEEPDISK))
    {
        // file doesn't exist, but is wanted
        DiskQTensor diskqt(this);  // should do everything in there
        // note - it's ok if the INMEM vs. DISK flags are incorrect in diskqt. They
        // aren't saved to the .dim file
    }
}

} // close namespace panache

