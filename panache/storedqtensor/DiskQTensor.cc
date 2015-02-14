/*! \file
 * \brief Three-index tensor storage on disk (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <cstdio> // for remove()
#include <sstream>

#include "panache/Exception.h"
#include "panache/Lapack.h"
#include "panache/ERI.h"
#include "panache/Flags.h"
#include "panache/storedqtensor/DiskQTensor.h"
#include "panache/storedqtensor/MemoryQTensor.h"
#include "panache/BasisSet.h"
#include "panache/FittingMetric.h"

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

void DiskQTensor::GenDFQso_(const SharedFittingMetric & fit,
                            const SharedBasisSet primary,
                            const SharedBasisSet auxiliary,
                            int nthreads)
{
    int maxpershell = primary->max_function_per_shell();
    int maxpershell2 = maxpershell*maxpershell;

    double * J = fit->get_metric();

    // default constructor = zero basis
    SharedBasisSet zero(new BasisSet);

    std::vector<SharedTwoBodyAOInt> eris;
    std::vector<const double *> eribuffers;
    std::vector<double *> A, B;

    int naux = StoredQTensor::naux();

    for(int i = 0; i < nthreads; i++)
    {
        eris.push_back(GetERI(auxiliary, zero, primary, primary));
        eribuffers.push_back(eris.back()->buffer());

        // temporary buffers
        A.push_back(new double[naux*maxpershell2]);
        B.push_back(new double[naux*maxpershell2]);
    }


    const int nprimshell = primary->nshell();
    const int nauxshell = auxiliary->nshell();

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
#endif
    for (int M = 0; M < nprimshell; M++)
    {
        int threadnum = 0;

#ifdef _OPENMP
        threadnum = omp_get_thread_num();
#endif

        int nm = primary->shell(M).nfunction();
        int mstart = primary->shell(M).function_index();
        int mend = mstart + nm;

        for (int N = 0; N <= M; N++)
        {
            int nn = primary->shell(N).nfunction();
            int nstart = primary->shell(N).function_index();
            //int nend = nstart + nn;

            for (int P = 0; P < nauxshell; P++)
            {
                int np = auxiliary->shell(P).nfunction();
                int pstart = auxiliary->shell(P).function_index();
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
            C_DGEMM('T','T',nm*nn, naux, naux, 1.0, B[threadnum], nm*nn, J, naux, 0.0,
                    A[threadnum], naux);


            // write to disk or store in memory
            if(N == M)
            {
                int iwrite = 1;
                for (int m0 = 0, m = mstart; m < mend; m0++, m++)
                    Write_(A[threadnum] + (m0*nm)*naux, iwrite++, calcindex(m, mstart));
            }
            else
            {
                for (int m0 = 0, m = mstart; m < mend; m0++, m++)
                    Write_(A[threadnum] + (m0*nn)*naux, nn, calcindex(m, nstart));
            }
        }
    }

    for(int i = 0; i < nthreads; i++)
    {
        delete [] A[i];
        delete [] B[i];
    }
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


void DiskQTensor::Finalize_(void)
{
}

void DiskQTensor::NoFinalize_(void)
{
}



} // close namespace panache

