#include "panache/Exception.h"
#include "panache/DFTensor2.h"

namespace panache
{

int DFTensor2::StoredQTensor::naux(void) const
{
    return naux_;
}

int DFTensor2::StoredQTensor::ndim1(void) const
{
    return ndim1_;
}

int DFTensor2::StoredQTensor::ndim2(void) const
{
    return ndim2_;
}

int DFTensor2::StoredQTensor::ndim12(void) const
{
    return ndim12_;
}

int DFTensor2::StoredQTensor::storesize(void) const
{
    return ndim12_*naux_;
}

int DFTensor2::StoredQTensor::packed(void) const
{
    return packed_;
}

int DFTensor2::StoredQTensor::byq(void) const
{
    return byq_;
}

int DFTensor2::StoredQTensor::calcindex(int i, int j) const
{
    if(!packed_)
        return (i*ndim2_+j);
    else if(i >= j)
        return ((i*(i+1))>>1) + j;
    else
        return ((j*(j+1))>>1) + i;
}

DFTensor2::StoredQTensor::StoredQTensor(int naux, int ndim1, int ndim2, bool packed, bool byq, QStorage storetype)
{
    naux_ = naux;
    ndim1_ = ndim1;
    ndim2_ = ndim2;
    storetype_ = storetype;
    packed_ = packed;
    byq_ = byq;

    if(packed && ndim1 != ndim2)
        throw RuntimeError("non square packed matrices?");

    ndim12_ = (packed ? (ndim1_ * (ndim2_+1))/2 : ndim1_*ndim2_);
}

DFTensor2::StoredQTensor::~StoredQTensor()
{
}

DFTensor2::QStorage DFTensor2::StoredQTensor::StoreType(void) const
{
    return storetype_;
}

void DFTensor2::StoredQTensor::Write(double * data, int i, int j)
{
    Write_(data, calcindex(i,j));
}

void DFTensor2::StoredQTensor::WriteByQ(double * data, int nq, int qstart)
{
    WriteByQ_(data, nq, qstart);
}

void DFTensor2::StoredQTensor::Read(double * data, int i, int j)
{
    Read_(data, calcindex(i,j));
}

int DFTensor2::StoredQTensor::ReadByQ(double * data, int nq, int qstart)
{
    if(qstart + nq >= naux_)
        nq = naux_-qstart;

    ReadByQ_(data, nq, qstart);
    return nq;
}

void DFTensor2::StoredQTensor::Reset(void)
{
    Reset_();
}

void DFTensor2::StoredQTensor::Clear(void)
{
    Clear_();
}

void DFTensor2::StoredQTensor::Init(void)
{
    Init_();
}



/*
    class DiskQTensor : public StoredQTensor
    {
    private:
        string filename_;
        std::unique_ptr<std::fstream> file_;

    protected:
        void Reset_(void)
        {
            if(file_)
            {
                file_->seekg(0);
                file_->seekp(0);
            }
        }

        void Write_(double * data, size_t nq, int ij)
        {
               // \todo write to file
        }

        void Read_(double * data, size_t nq, int ij)
        {
               // \todo write to file
        }

    public:
        DiskQTensor(int naux, int ndim1, int ndim2, bool packed, const string & filename)
            : StoredQTensor(naux, ndim1, ndim2, packed, QStorage::ONDISK)
        {
            filename_ = filename;
        }

        void OpenFile(void)
        {
            if(file_is_open())
                return;


            if(filename.length() == 0)
                throw RuntimeError("Error - no file specified!");

            file_ = std::unique_ptr<std::fstream>(new std::fstream(filename_.c_str()),
                                                  std::fstream::in | std::fstream::out |
                                                  std::fstream::binary | std::fstream::trunc );
            if(!file_->is_open())
                throw RuntimeError(filename_);

            file_->exceptions(std::fstream::failbit | std::fstream::badbit | std::fstream::eofbit);
            curij_ = 0;
        }

        void CloseFile(void)
        {
            if(file_ && file_->is_open())
            {
                file_->close();
                matfile_.reset();
            }
        }


    };
*/


//////////////////////////////
// MemoryQTensor
//////////////////////////////
void DFTensor2::MemoryQTensor::Reset_(void)
{
    // nothing needed
}

void DFTensor2::MemoryQTensor::Write_(double * data, int ij)
{
    if(byq())
    {
        for(int q = 0; q < naux(); q++)
        {
            data_[q*ndim12()+ij] = data[q];
        }
    }
    else
    {
        double * start = data_.get() + ij * naux();
        std::copy(data, data+naux(), start);
    }
}

void DFTensor2::MemoryQTensor::WriteByQ_(double * data, int nq, int qstart)
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

void DFTensor2::MemoryQTensor::Read_(double * data, int ij)
{
    // index ij is given by calling function and takes into account packing
    if(byq())
    {
        for(int q0 = 0; q0 < naux(); q0++)
            data[q0]
                = data_[q0*ndim12()+ij];
    }
    else
    {
        double * start = data_.get() + ij * naux();
        std::copy(start, start+naux(), data);
    }
}

void DFTensor2::MemoryQTensor::ReadByQ_(double * data, int nq, int qstart)
{
    if(packed())
    {
        if(byq())
        {
            for(int q0 = 0, q = qstart; q0 < nq; q++, q0++)
                for(int m = 0; m < ndim1(); m++)
                    for(int n = 0; n <= m; n++)
                        data[q0*ndim1()*ndim2() + m*ndim2() + n]
                            = data[q0*ndim1()*ndim2() + n*ndim1() + m]
                              = data_[q*ndim12()+calcindex(m,n)];
        }
        else
        {
            for(int q0 = 0, q = qstart; q0 < nq; q++, q0++)
                for(int m = 0; m < ndim1(); m++)
                    for(int n = 0; n <= m; n++)
                        data[q0*ndim1()*ndim2() + m*ndim2() + n]
                            = data[q0*ndim1()*ndim2() + n*ndim1() + m]
                              = data_[calcindex(m,n)*naux() + q];
        }
    }
    else
    {
        if(byq())
        {
            for(int q0 = 0, q = qstart; q0 < nq; q++, q0++)
                for(int m = 0; m < ndim1(); m++)
                    for(int n = 0; n < ndim2(); n++)
                        data[q0*ndim1()*ndim2() + m*ndim2() + n]
                            = data_[q*ndim12() + calcindex(m,n)];
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
}

void DFTensor2::MemoryQTensor::Clear_(void)
{
    data_.reset();
}

void DFTensor2::MemoryQTensor::Init_(void)
{
    if(!data_)
        data_ = std::unique_ptr<double []>(new double[storesize()]);
}

DFTensor2::MemoryQTensor::MemoryQTensor(int naux, int ndim1, int ndim2, bool packed, bool byq)
    : StoredQTensor(naux, ndim1, ndim2, packed, byq, QStorage::INMEM)
{
}




std::unique_ptr<DFTensor2::StoredQTensor> DFTensor2::StoredQTensorFactory(int naux, int ndim1, int ndim2, bool packed, bool byq, QStorage storetype)
{
    if(storetype == QStorage::INMEM)
        return std::unique_ptr<DFTensor2::StoredQTensor>(new MemoryQTensor(naux, ndim1, ndim2, packed, byq));
    else
        throw RuntimeError("No StoredQTensor for that type");
}

} // close namespace panache

