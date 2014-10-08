/*! \file
 * \brief Three-index tensor storage on disk (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/Exception.h"
#include "panache/storedqtensor/DiskQTensor.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace panache
{

void DiskQTensor::OpenFile_(void)
{
    if(file_ && file_->is_open())
        return;

    if(name().length() == 0)
        throw RuntimeError("Error - no file specified!");

    file_ = std::unique_ptr<std::fstream>(new std::fstream(name().c_str(),
                                          std::fstream::in | std::fstream::out |
                                          std::fstream::binary | std::fstream::trunc ));
    if(!file_->is_open())
        throw RuntimeError(name());

    file_->exceptions(std::fstream::failbit | std::fstream::badbit | std::fstream::eofbit);
}

void DiskQTensor::CloseFile_(void)
{
    if(file_ && file_->is_open())
    {
        file_->close();
        file_.reset();
    }
}

void DiskQTensor::Reset_(void)
{
    file_->seekg(0);
    file_->seekp(0);
}

void DiskQTensor::Write_(double * data, int nij, int ijstart)
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

void DiskQTensor::WriteByQ_(double * data, int nq, int qstart)
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

void DiskQTensor::Read_(double * data, int nij, int ijstart)
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

void DiskQTensor::ReadByQ_(double * data, int nq, int qstart)
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
                file_->read(reinterpret_cast<char *>(data+q0*ndim12() + calcindex(m,n)), sizeof(double));
            }
        }
    }
}

void DiskQTensor::Clear_(void)
{
    //! \todo Erase file
    CloseFile_();
}

void DiskQTensor::Init_(void)
{
    OpenFile_();
}


DiskQTensor::DiskQTensor()
{
}


} // close namespace panache

