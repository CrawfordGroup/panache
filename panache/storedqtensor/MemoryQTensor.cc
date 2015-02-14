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
#include <iostream>

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


void MemoryQTensor::WriteByQ_(double * data, int nq, int qstart)
{
    int indim12 = ndim12();

    if(byq())
    {
        std::copy(data,
                  data + nq*indim12,
                  data_.get()+qstart*indim12);
    }
    else
    {
        int inaux = naux();

        for(int q = 0, qoff = 0; q < nq; q++, qoff += indim12)
            for(int ij = 0, ijoff = 0; ij < indim12; ij++, ijoff += inaux)
                data_[ijoff+qstart+q] = data[qoff+ij];
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

void MemoryQTensor::GenDFQso_(const SharedFittingMetric & fit,
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

void MemoryQTensor::Finalize_(void)
{
}

void MemoryQTensor::NoFinalize_(void)
{
}

} // close namespace panache

