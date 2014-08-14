#include "panache/Exception.h"
#include "panache/DFTensor.h"
namespace panache
{


/////////////////////////////
// STOREDQTENSOR BASE CLASS
/////////////////////////////
int DFTensor::StoredQTensor::naux(void) const
{
    return naux_;
}

int DFTensor::StoredQTensor::ndim1(void) const
{
    return ndim1_;
}

int DFTensor::StoredQTensor::ndim2(void) const
{
    return ndim2_;
}

int DFTensor::StoredQTensor::ndim12(void) const
{
    return ndim12_;
}

int DFTensor::StoredQTensor::storesize(void) const
{
    return ndim12_*naux_;
}

int DFTensor::StoredQTensor::packed(void) const
{
    return (storeflags_ & QSTORAGE_PACKED);
}

int DFTensor::StoredQTensor::byq(void) const
{
    return (storeflags_ & QSTORAGE_BYQ);
}

int DFTensor::StoredQTensor::calcindex(int i, int j) const
{
    if(!packed())
        return (i*ndim2_+j);
    else if(i >= j)
        return ((i*(i+1))>>1) + j;
    else
        return ((j*(j+1))>>1) + i;
}

DFTensor::StoredQTensor::StoredQTensor(int naux, int ndim1, int ndim2, int storeflags)
{
    naux_ = naux;
    ndim1_ = ndim1;
    ndim2_ = ndim2;
    storeflags_ = storeflags;

    if(packed() && ndim1 != ndim2)
        throw RuntimeError("non square packed matrices?");

    ndim12_ = (packed() ? (ndim1_ * (ndim2_+1))/2 : ndim1_*ndim2_);
}

DFTensor::StoredQTensor::~StoredQTensor()
{
}

int DFTensor::StoredQTensor::StoreFlags(void) const
{
    return storeflags_;
}

void DFTensor::StoredQTensor::Write(double * data, int nij, int ijstart)
{
    Write_(data, nij, ijstart);
}

void DFTensor::StoredQTensor::WriteByQ(double * data, int nq, int qstart)
{
    WriteByQ_(data, nq, qstart);
}

int DFTensor::StoredQTensor::Read(double * data, int nij, int ijstart)
{
    if(ijstart + nij >= ndim12())
        nij = ndim12() - ijstart;

    Read_(data, nij, ijstart);
    return nij;
}

int DFTensor::StoredQTensor::ReadByQ(double * data, int nq, int qstart)
{
    if(qstart + nq >= naux_)
        nq = naux_-qstart;

    ReadByQ_(data, nq, qstart);
    return nq;
}

void DFTensor::StoredQTensor::Reset(void)
{
    Reset_();
}

void DFTensor::StoredQTensor::Clear(void)
{
    Clear_();
}

void DFTensor::StoredQTensor::Init(void)
{
    Init_();
}


CumulativeTime & DFTensor::StoredQTensor::GenTimer(void)
{
    return gen_timer_; 
}

CumulativeTime & DFTensor::StoredQTensor::GetQBatchTimer(void)
{
    return getqbatch_timer_;
}

CumulativeTime & DFTensor::StoredQTensor::GetBatchTimer(void)
{
    return getijbatch_timer_;
}


//////////////////////////////
// DiskQTensor
//////////////////////////////
void DFTensor::DiskQTensor::OpenFile_(void)
{
    if(file_ && file_->is_open())
        return;

    if(filename_.length() == 0)
        throw RuntimeError("Error - no file specified!");

    file_ = std::unique_ptr<std::fstream>(new std::fstream(filename_.c_str(),
                                          std::fstream::in | std::fstream::out |
                                          std::fstream::binary | std::fstream::trunc ));
    if(!file_->is_open())
        throw RuntimeError(filename_);

    file_->exceptions(std::fstream::failbit | std::fstream::badbit | std::fstream::eofbit);
}

void DFTensor::DiskQTensor::CloseFile_(void)
{
    if(file_ && file_->is_open())
    {
        file_->close();
        file_.reset();
    }
}

void DFTensor::DiskQTensor::Reset_(void)
{
    file_->seekg(0);
    file_->seekp(0);
}

void DFTensor::DiskQTensor::Write_(double * data, int nij, int ijstart)
{
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
        if(byq())
        {
            for(int q = 0; q < naux(); q++)
            {
                file_->seekp(sizeof(double)*(q*ndim12()+ijstart), std::ios_base::beg);

                for(int ij = 0; ij < nij; ij++)
                    file_->write(reinterpret_cast<const char *>(data+ij*naux()+q), sizeof(double));
            }
        }
        else
        {
            file_->seekp(sizeof(double)*(ijstart * naux()), std::ios_base::beg);
            file_->write(reinterpret_cast<const char *>(data), nij*naux()*sizeof(double));
        }
    }
}

void DFTensor::DiskQTensor::WriteByQ_(double * data, int nq, int qstart)
{
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
        if(byq())
        {
            file_->seekp(sizeof(double)*(qstart*ndim12()), std::ios_base::beg);
            file_->write(reinterpret_cast<const char *>(data), nq*ndim12()*sizeof(double));
        }
        else
        {
            for(int q = 0; q < nq; q++)
            {
                for(int ij = 0; ij < ndim12(); ij++)
                {
                    file_->seekp(sizeof(double)*(ij*naux()+qstart+q), std::ios_base::beg);
                    file_->write(reinterpret_cast<const char *>(data+q*ndim12()+ij), sizeof(double));
                }
            }
        }
    }
}

void DFTensor::DiskQTensor::Read_(double * data, int nij, int ijstart)
{
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
        if(byq())
        {
            for(int q = 0; q < naux(); q++)
            {
                file_->seekg(sizeof(double)*(q*ndim12()+ijstart), std::ios_base::beg);

                for(int ij = 0; ij < nij; ij++)
                    file_->read(reinterpret_cast<char *>(data+ij*naux()+q), sizeof(double));
            }
        }
        else
        {
            file_->seekg(sizeof(double)*(ijstart*naux()), std::ios_base::beg);
            file_->read(reinterpret_cast<char *>(data), sizeof(double)*nij*naux());
        }
    }
}

void DFTensor::DiskQTensor::ReadByQ_(double * data, int nq, int qstart)
{
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
        if(byq())
        {
            file_->seekg(sizeof(double)*(qstart*ndim12()), std::ios_base::beg);
            file_->read(reinterpret_cast<char *>(data), sizeof(double)*nq*ndim12());
        }
        else
        {
            for(int q0 = 0, q = qstart; q0 < nq; q++, q0++)
            for(int m = 0; m < ndim1(); m++)
            for(int n = 0; n < ndim2(); n++)
            {
                file_->seekg(sizeof(double)*(calcindex(m,n)*naux()+q), std::ios_base::beg);
                file_->read(reinterpret_cast<char *>(data+q0*ndim1()*ndim2() + m*ndim2()+n), sizeof(double));
            }
        }
    }
}

void DFTensor::DiskQTensor::Clear_(void)
{
    //! \todo Erase file
    CloseFile_();
}

void DFTensor::DiskQTensor::Init_(void)
{
    OpenFile_();
}


DFTensor::DiskQTensor::DiskQTensor(int naux, int ndim1, int ndim2, int storeflags, const std::string & filename)
            : StoredQTensor(naux, ndim1, ndim2, storeflags)
{
    filename_ = filename;
}



//////////////////////////////
// MemoryQTensor
//////////////////////////////
void DFTensor::MemoryQTensor::Reset_(void)
{
    // nothing needed
}

void DFTensor::MemoryQTensor::Write_(double * data, int nij, int ijstart)
{
    if(byq())
    {
        for(int q = 0; q < naux(); q++)
        for(int ij = 0; ij < nij; ij++)
            data_[q*ndim12()+ijstart+ij] = data[ij*naux()+q];
    }
    else
    {
        double * start = data_.get() + ijstart * naux();
        std::copy(data, data+naux(), start);
    }
}

void DFTensor::MemoryQTensor::WriteByQ_(double * data, int nq, int qstart)
{
    if(byq())
    {
        std::copy(data,
                  data + nq*ndim12(),
                  data_.get()+qstart*ndim12());
    }
    else
    {
        for(int q = 0; q < nq; q++)
        {
            for(int ij = 0; ij < ndim12(); ij++)
                data_[ij*naux()+qstart+q] = data[q*ndim12()+ij];
        }
    }
}

void DFTensor::MemoryQTensor::Read_(double * data, int nij, int ijstart)
{
    // index ij is given by calling function and takes into account packing
    if(byq())
    {
        for(int q = 0; q < naux(); q++)
        for(int ij = 0; ij < nij; ij++)
            data[ij*naux()+q] = data_[q*ndim12()+ijstart+ij];
    }
    else
    {
        double * start = data_.get() + ijstart * naux();
        std::copy(start, start+nij*naux(), data);
    }
}

void DFTensor::MemoryQTensor::ReadByQ_(double * data, int nq, int qstart)
{
    if(byq())
    {
        double * start = data_.get()+qstart*ndim12();
        std::copy(start, start+nq*ndim12(), data);
    }
    else
    {
        for(int q0 = 0, q = qstart; q0 < nq; q++, q0++)
        for(int m = 0; m < ndim1(); m++)
        for(int n = 0; n < ndim2(); n++)
            data[q0*ndim1()*ndim2() + m*ndim2() + n]
                = data_[calcindex(m,n)*naux() + q];
    }
}

void DFTensor::MemoryQTensor::Clear_(void)
{
    data_.reset();
}

void DFTensor::MemoryQTensor::Init_(void)
{
    if(!data_)
        data_ = std::unique_ptr<double []>(new double[storesize()]);
}

DFTensor::MemoryQTensor::MemoryQTensor(int naux, int ndim1, int ndim2, int storeflags)
    : StoredQTensor(naux, ndim1, ndim2, storeflags)
{
}




std::unique_ptr<DFTensor::StoredQTensor> DFTensor::StoredQTensorFactory(int naux, int ndim1, int ndim2, 
                                                                        int storeflags, const std::string & name)
{
    if(name == "")
        throw RuntimeError("NO NAME SPECIFIED");

    if(storeflags & QSTORAGE_ONDISK)
    {
        std::string filename(directory_);
        filename.append("/");
        filename.append(name);
        return std::unique_ptr<DFTensor::StoredQTensor>(new DiskQTensor(naux, ndim1, ndim2, storeflags, filename));
    }
    else
        return std::unique_ptr<DFTensor::StoredQTensor>(new MemoryQTensor(naux, ndim1, ndim2, storeflags));
}

} // close namespace panache

