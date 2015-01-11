/*! \file
 * \brief Three-index tensor storage on disk (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <cstdio> // for remove()
#include <sstream>

#include "panache/Exception.h"
#include "panache/Flags.h"
#include "panache/storedqtensor/DiskQTensor.h"
#include "panache/storedqtensor/MemoryQTensor.h"

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


void DiskQTensor::OpenForReadWrite_(void)
{
    file_ = std::unique_ptr<std::fstream>(new std::fstream(filename_.c_str(),
                                          std::fstream::out | std::fstream::in |
                                          std::fstream::trunc | std::fstream::binary));

    if(!file_->is_open())
        throw RuntimeError(std::string("Unable to open file ") + filename_);

    file_->exceptions(std::fstream::failbit | std::fstream::badbit | std::fstream::eofbit);
}


bool DiskQTensor::OpenForRead_(bool required)
{
    file_ = std::unique_ptr<std::fstream>(new std::fstream(filename_.c_str(),
                                          std::fstream::in | std::fstream::binary));

    if(!file_->is_open())
    {
        if(required)
            throw RuntimeError(std::string("Unable to open file ") + filename_);
        else
        {
            file_.reset();
            return false;
        }
    }

    file_->exceptions(std::fstream::failbit | std::fstream::badbit | std::fstream::eofbit);
    return true;
}


void DiskQTensor::Init_(void)
{
    WriteDimFile_();
    // Nothing else to really do here. Reserve a big file?
}


DiskQTensor::DiskQTensor(int storeflags, const std::string & name, const std::string & directory) 
             : LocalQTensor(storeflags, name, directory)
{
    if(file_ && file_->is_open())
        return;

    if(filename_.length() == 0 || name.length() == 0)
        throw RuntimeError("Error - no file specified!");


    if(existed_ && (storeflags & QSTORAGE_READDISK))
    {
        OpenForRead_(true);  // false = not required
        ReadDimFile_();
        markfilled();
    }

    if(!(storeflags & QSTORAGE_READDISK) || !existed_)
        OpenForReadWrite_();
}


void DiskQTensor::ReadDimFile_(void)
{
    std::ifstream dim(dimfilename_.c_str());

    if(!dim.is_open())
        throw RuntimeError(std::string("Unable to open file ") + dimfilename_);

    dim.exceptions(std::fstream::failbit | std::fstream::badbit | std::fstream::eofbit);

    dim >> f_naux_ >> f_ndim1_ >> f_ndim2_ >> f_ndim12_ >> f_packed_ >> f_byq_;

    std::stringstream ss;
 
    // careful. byq() and ibyq are both ints and would represent QSTORAGE_BYQ, etc, not just a simple bool
    if(f_byq_ != byq())
    {
        ss << "Tensor " << name() << " does not match orientation from file " << dimfilename_
           << " Here: " << byq() << " disk: " << f_byq_ << "\n";
        throw RuntimeError(ss.str());
    }

    // same here
    if(f_packed_ != packed())
    {
        ss << "Tensor " << name() << " does not match packed-ness from file " << dimfilename_
           << " Here: " << packed() << " disk: " << f_packed_ << "\n";
        throw RuntimeError(ss.str());
    }

    Init(f_naux_, f_ndim1_, f_ndim2_);
}

void DiskQTensor::WriteDimFile_(void)
{
    std::ofstream dim(dimfilename_.c_str(), std::ofstream::trunc);

    if(!dim.is_open())
        throw RuntimeError(std::string("Unable to open file ") + dimfilename_);

    dim.exceptions(std::fstream::failbit | std::fstream::badbit | std::fstream::eofbit);

    // careful. byq() and packed() are both ints and would represent QSTORAGE_BYQ, etc, not just a simple bool
    dim << naux() << " " << ndim1() << " " << ndim2() << " "
        << ndim12() << " " << packed() << " "
        << byq() << " END";  //Sorry, the "END" is a cheap hack so that ReadDimFile_ doesn't
                             // throw with EOF
}


DiskQTensor::~DiskQTensor()
{
    if(file_ && file_->is_open())
    {
        file_->close();
        file_.reset();
    }

    // Erase the file
    if(!(storeflags() & QSTORAGE_KEEPDISK))
    {
        std::remove(filename_.c_str());
        std::remove(dimfilename_.c_str());
    }
}


DiskQTensor::DiskQTensor(MemoryQTensor * memqt) 
                 : DiskQTensor(memqt->storeflags(), memqt->name(), memqt->directory())
{
    // initialize sizes
    // from StoredQTensor base class
    Init(*memqt);

    int inaux = naux();
    int indim12 = ndim12();

    // do in blocks
    if(byq())
    {
        std::unique_ptr<double[]> buf(new double[indim12]);
        double * bufptr = buf.get();
  
        for(int i = 0; i < inaux; i++)
        {
            memqt->ReadByQ(bufptr, 1, i);
            WriteByQ_(bufptr, 1, i);
        }
    }
    else
    {
        std::unique_ptr<double[]> buf(new double[inaux]);
        double * bufptr = buf.get();

        for(int i = 0; i < indim12; i++)
        {
            memqt->Read(bufptr, 1, i);
            Write_(bufptr, 1, i);
        }
    }
}





} // close namespace panache

