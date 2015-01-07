/*! \file
 * \brief Three-index tensor storage on disk (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <cstdio> // for remove()

#include "panache/Exception.h"
#include "panache/storedqtensor/DiskQTensor.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace panache
{

void DiskQTensor::Write_(double * data, int nij, int ijstart)
{
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
        int inaux = naux();

        if(byq())
        {
            int indim12 = ndim12();

            for(int q = 0, qoff = 0; q < inaux; q++, qoff += indim12)
            {
                file_->seekp(sizeof(double)*(qoff+ijstart), std::ios_base::beg);

                for(int ij = 0, ijoff = 0; ij < nij; ij++, ijoff += inaux)
                    file_->write(reinterpret_cast<const char *>(data+ijoff+q), sizeof(double));
            }
        }
        else
        {
            file_->seekp(sizeof(double)*(ijstart * inaux), std::ios_base::beg);
            file_->write(reinterpret_cast<const char *>(data), nij*inaux*sizeof(double));
        }
    }
}

void DiskQTensor::WriteByQ_(double * data, int nq, int qstart)
{
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
        int indim12 = ndim12();

        if(byq())
        {
            file_->seekp(sizeof(double)*(qstart*indim12), std::ios_base::beg);
            file_->write(reinterpret_cast<const char *>(data), nq*indim12*sizeof(double));
        }
        else
        {
            int inaux = naux();

            for(int q = 0, qoff = 0; q < nq; q++, qoff += indim12)
            {
                for(int ij = 0, ijoff = 0; ij < indim12; ij++, ijoff += inaux)
                {
                    file_->seekp(sizeof(double)*(ijoff+qstart+q), std::ios_base::beg);
                    file_->write(reinterpret_cast<const char *>(data+qoff+ij), sizeof(double));
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
        int inaux = naux();

        if(byq())
        {
            int indim12 = ndim12();

            for(int q = 0, qoff = 0; q < inaux; q++, qoff += indim12)
            {
                file_->seekg(sizeof(double)*(qoff+ijstart), std::ios_base::beg);

                for(int ij = 0, ijoff = 0; ij < nij; ij++, ijoff += inaux)
                    file_->read(reinterpret_cast<char *>(data+ijoff+q), sizeof(double));
            }
        }
        else
        {
            file_->seekg(sizeof(double)*(ijstart*inaux), std::ios_base::beg);
            file_->read(reinterpret_cast<char *>(data), sizeof(double)*nij*inaux);
        }
    }
}

void DiskQTensor::ReadByQ_(double * data, int nq, int qstart)
{
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
        int indim12 = ndim12();

        if(byq())
        {
            file_->seekg(sizeof(double)*(qstart*indim12), std::ios_base::beg);
            file_->read(reinterpret_cast<char *>(data), sizeof(double)*nq*indim12);
        }
        else
        {
            int inaux = naux();

            for(int q = 0, qoff = 0; q < nq; q++, qoff += indim12)
            for(int ij = 0, ijoff = 0; ij < indim12; ij++, ijoff += inaux)
            {
                file_->seekg(sizeof(double)*(ijoff+qstart+q), std::ios_base::beg);
                file_->read(reinterpret_cast<char *>(data+qoff+ij), sizeof(double));
            }
        }
    }
}

void DiskQTensor::Init_(void)
{
    filename_ = directory_;
    filename_.append("/");
    filename_.append(name());

    if(file_ && file_->is_open())
        return;

    if(filename_.length() == 0 || name().length() == 0)
        throw RuntimeError("Error - no file specified!");

    file_ = std::unique_ptr<std::fstream>(new std::fstream(filename_.c_str(),
                                          std::fstream::in | std::fstream::out |
                                          std::fstream::binary | std::fstream::trunc ));
    if(!file_->is_open())
        throw RuntimeError(std::string("Unable to open file ") + filename_);

    file_->exceptions(std::fstream::failbit | std::fstream::badbit | std::fstream::eofbit);
}


DiskQTensor::DiskQTensor(const std::string & directory) : directory_(directory)
{
}


DiskQTensor::~DiskQTensor()
{
    if(file_ && file_->is_open())
    {
        file_->close();
        file_.reset();
    }

    // erase the file
    std::remove(filename_.c_str());
}



} // close namespace panache

