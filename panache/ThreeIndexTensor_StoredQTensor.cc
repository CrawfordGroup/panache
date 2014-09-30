#include <fstream>

#include "panache/Exception.h"
#include "panache/ThreeIndexTensor.h"
#include "panache/Lapack.h"
#include "panache/BasisSet.h"
#include "panache/FittingMetric.h"
#include "panache/Timing.h"
#include "panache/Lapack.h"
#include "panache/ERI.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace panache
{


/////////////////////////////
// STOREDQTENSOR BASE CLASS
/////////////////////////////
int ThreeIndexTensor::StoredQTensor::naux(void) const
{
    return naux_;
}

int ThreeIndexTensor::StoredQTensor::ndim1(void) const
{
    return ndim1_;
}

int ThreeIndexTensor::StoredQTensor::ndim2(void) const
{
    return ndim2_;
}

int ThreeIndexTensor::StoredQTensor::ndim12(void) const
{
    return ndim12_;
}

int ThreeIndexTensor::StoredQTensor::storesize(void) const
{
    return ndim12_*naux_;
}

int ThreeIndexTensor::StoredQTensor::packed(void) const
{
    return (storeflags_ & QSTORAGE_PACKED);
}

int ThreeIndexTensor::StoredQTensor::byq(void) const
{
    return (storeflags_ & QSTORAGE_BYQ);
}

const std::string & ThreeIndexTensor::StoredQTensor::name(void) const
{
    return name_;
}
int ThreeIndexTensor::StoredQTensor::calcindex(int i, int j) const
{
    if(!packed())
        return (i*ndim2_+j);
    else if(i >= j)
        return ((i*(i+1))>>1) + j;
    else
        return ((j*(j+1))>>1) + i;
}

void ThreeIndexTensor::StoredQTensor::Init(int naux, int ndim1, int ndim2, int storeflags, const std::string & name)
{
    naux_ = naux;
    ndim1_ = ndim1;
    ndim2_ = ndim2;
    storeflags_ = storeflags;

    if(packed() && ndim1 != ndim2)
        throw RuntimeError("non square packed matrices?");

    ndim12_ = (packed() ? (ndim1_ * (ndim2_+1))/2 : ndim1_*ndim2_);
    name_ = name;

    Init_();
}

ThreeIndexTensor::StoredQTensor::~StoredQTensor()
{
}

int ThreeIndexTensor::StoredQTensor::StoreFlags(void) const
{
    return storeflags_;
}

int ThreeIndexTensor::StoredQTensor::Read(double * data, int nij, int ijstart)
{
    if(ijstart + nij >= ndim12())
        nij = ndim12() - ijstart;

    Read_(data, nij, ijstart);
    return nij;
}

int ThreeIndexTensor::StoredQTensor::ReadByQ(double * data, int nq, int qstart)
{
    if(qstart + nq >= naux_)
        nq = naux_-qstart;

    ReadByQ_(data, nq, qstart);
    return nq;
}

void ThreeIndexTensor::StoredQTensor::Reset(void)
{
    Reset_();
}

void ThreeIndexTensor::StoredQTensor::Clear(void)
{
    Clear_();
}

ThreeIndexTensor::StoredQTensor::StoredQTensor(void)
{
}


CumulativeTime & ThreeIndexTensor::StoredQTensor::GenTimer(void)
{
    return gen_timer_;
}

CumulativeTime & ThreeIndexTensor::StoredQTensor::GetQBatchTimer(void)
{
    return getqbatch_timer_;
}

CumulativeTime & ThreeIndexTensor::StoredQTensor::GetBatchTimer(void)
{
    return getijbatch_timer_;
}

void ThreeIndexTensor::StoredQTensor::GenDFQso(const std::shared_ptr<FittingMetric> & fit,
                                     const SharedBasisSet primary,
                                     const SharedBasisSet auxiliary,
                                     int nthreads)
{
    GenDFQso_(fit, primary, auxiliary, nthreads);
}

void ThreeIndexTensor::StoredQTensor::GenCHQso(const SharedBasisSet primary,
                                     double delta,
                                     int storeflags,
                                     int nthreads)
{
    GenCHQso_(primary, delta, storeflags, nthreads);
}


void ThreeIndexTensor::StoredQTensor::Transform(const std::vector<TransformMat> & left,
                                        const std::vector<TransformMat> & right,
                                        std::vector<StoredQTensor *> results,
                                        int nthreads)
{
    Transform_(left, right, results, nthreads);
}




//////////////////////////////
// LocalQTensor
//////////////////////////////


ThreeIndexTensor::LocalQTensor::LocalQTensor()
{
}



void ThreeIndexTensor::LocalQTensor::GenDFQso_(const std::shared_ptr<FittingMetric> & fit,
                                     const SharedBasisSet primary,
                                     const SharedBasisSet auxiliary,
                                     int nthreads)
{
#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
#endif

    int maxpershell = primary->max_function_per_shell();
    int maxpershell2 = maxpershell*maxpershell;

    double * J = fit->get_metric();

    // default constructor = zero basis
    SharedBasisSet zero(new BasisSet);

    std::vector<std::shared_ptr<TwoBodyAOInt>> eris;
    std::vector<const double *> eribuffers;
    std::vector<double *> A, B;

    int naux = StoredQTensor::naux();

    for(int i = 0; i < nthreads; i++)
    {
        eris.push_back(GetERI(auxiliary, zero, primary, primary));
        eribuffers.push_back(eris[eris.size()-1]->buffer());

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

#ifdef PANACHE_TIMING
    tim.Stop();
    GenTimer().AddTime(tim);
#endif

}


void ThreeIndexTensor::LocalQTensor::ComputeDiagonal_(TwoBodyAOInt * integral, double * target)
{
    const double* buffer = integral->buffer();
    SharedBasisSet basis = integral->basis();

    for (int M = 0; M < basis->nshell(); M++) 
    {
        int nM = basis->shell(M).nfunction();
        int mstart = basis->shell(M).function_index();

        for (int N = 0; N < basis->nshell(); N++) 
        {
            int nN = basis->shell(N).nfunction();
            int nstart = basis->shell(N).function_index();

            integral->compute_shell(M,N,M,N);

            for (int om = 0; om < nM; om++) {
                for (int on = 0; on < nN; on++) {
                    target[(om + mstart) * basis->nbf() + (on + nstart)] =
                        buffer[om * nN * nM * nN + on * nM * nN + om * nN + on];
                }
            }
        }
    }
}


void ThreeIndexTensor::LocalQTensor::ComputeRow_(TwoBodyAOInt * integral, int row, double* target)
{
    const double* buffer = integral->buffer();
    SharedBasisSet basis = integral->basis();

    int r = row / basis->nbf();
    int s = row % basis->nbf();
    int R = basis->function_to_shell(r);
    int S = basis->function_to_shell(s);

    int nR = basis->shell(R).nfunction();
    int nS = basis->shell(S).nfunction();
    int rstart = basis->shell(R).function_index();
    int sstart = basis->shell(S).function_index();

    int oR = r - rstart;
    int os = s - sstart;

    for (int M = 0; M < basis->nshell(); M++)
    {
        int nM = basis->shell(M).nfunction();
        int mstart = basis->shell(M).function_index();

        for (int N = 0; N < basis->nshell(); N++)
        {
            integral->compute_shell(M,N,R,S);

            int nN = basis->shell(N).nfunction();
            int nstart = basis->shell(N).function_index();

            for (int om = 0; om < nM; om++) {
                for (int on = 0; on < nN; on++) {
                    target[(om + mstart) * basis->nbf() + (on + nstart)] =
                        buffer[om * nN * nR * nS + on * nR * nS + oR * nS + os];
                }
            }
        }
    }
}



void ThreeIndexTensor::LocalQTensor::GenCHQso_(const SharedBasisSet primary,
                                               double delta,
                                               int storeflags,
                                               int nthreads)
{
#ifdef PANACHE_TIMING
    Timer tim;
    tim.Start();
#endif

    std::vector<std::shared_ptr<TwoBodyAOInt>> eris;

    for(int i = 0; i < nthreads; i++)
        eris.push_back(GetERI(primary, primary, primary, primary));

    int nQ = 0;
    int n = primary->nbf();
    int n2 = n*n;

    double * diag = new double[n2];
    ComputeDiagonal_(eris[0].get(), diag);


    // Temporary cholesky factor
    std::vector<double*> L;

    // List of selected pivots
    std::vector<int> pivots;
 
    while(nQ < n2)
    {
        int pivot = 0;
        double Dmax = diag[0];
        for(int P = 0; P < n2; P++)
        {
            if(Dmax < diag[P])
            {
                Dmax = diag[P];
                pivot = P;
            }
        }

        if(Dmax < delta || Dmax < 0.0) break;

        pivots.push_back(pivot);
        double L_QQ = sqrt(Dmax);

        L.push_back(new double[n2]);
        ComputeRow_(eris[0].get(), pivot, L[nQ]);

        // [(m|Q) - L_m^P L_Q^P]
        for (int P = 0; P < nQ; P++) {
            C_DAXPY(n2,-L[P][pivots[nQ]],L[P],1,L[nQ],1);
        }

        // 1/L_QQ [(m|Q) - L_m^P L_Q^P]
        C_DSCAL(n2, 1.0 / L_QQ, L[nQ], 1);

        // Zero the upper triangle
        for (size_t P = 0; P < pivots.size(); P++) {
            L[nQ][pivots[P]] = 0.0;
        }

        // Set the pivot factor
        L[nQ][pivot] = L_QQ;

        // Update the Schur complement diagonal
        for (int P = 0; P < n2; P++) {
            diag[P] -= L[nQ][P] * L[nQ][P];
        }

        // Force truly zero elements to zero
        for (size_t P = 0; P < pivots.size(); P++) {
            diag[pivots[P]] = 0.0;
        }

        nQ++;
    }

    // copy to memory now that we have the sizes
    StoredQTensor::Init(nQ, n, n, storeflags | QSTORAGE_BYQ, "qso");

    for(int i = 0; i < nQ; i++)
    {
        WriteByQ_(L[i], 1, i);
        delete [] L[i];
    }
     
}


void ThreeIndexTensor::LocalQTensor::Transform_(const std::vector<TransformMat> & left,
                                        const std::vector<TransformMat> & right,
                                        std::vector<StoredQTensor *> results,
                                        int nthreads)
{
    int naux = StoredQTensor::naux();
    int ndim1 = StoredQTensor::ndim1();
    int ndim2 = StoredQTensor::ndim2();
    int ndim12 = StoredQTensor::ndim12();

    if(left.size() != right.size())
        throw RuntimeError("Error - imbalanced left & right transformation matrix sizes!");
    if(left.size() != results.size())
        throw RuntimeError("Error - not enough results in vector");

    // find the max dimension of the transformation matrices
    int maxl = 0;
    int maxr = 0;
    for(const auto & it : left)
        maxl = (it.second > maxl) ? it.second : maxl;
    for(const auto & it : right)
        maxr = (it.second > maxr) ? it.second : maxr;

    // temporary space
    double * qe = new double[ndim1*ndim2*nthreads];   // expanded q
    double * qc = new double[maxl*ndim2*nthreads];     // C(t) Q
    double * cqc = new double[maxl*maxr*nthreads];    // C(t) Q C

    double * qp = qe;
    if(packed())
        qp = new double[ndim12*nthreads];         // packed q


    #ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads)
    #endif
    for(int q = 0; q < naux; q++)
    {
        int threadnum = 0;

        #ifdef _OPENMP
            threadnum = omp_get_thread_num();
        #endif

        for(size_t i = 0; i < left.size(); i++)
        {

            #ifdef PANACHE_TIMING
                Timer tim;
                tim.Start();
            #endif

            int lncols = left[i].second;
            int rncols = right[i].second;
            double * lptr = left[i].first;
            double * rptr = right[i].first;

            LocalQTensor * qout = dynamic_cast<LocalQTensor *>(results[i]);
            if(qout == nullptr)
                throw RuntimeError("Cannot transform LocalQTensor into another type!");

            double * myqe = qe + threadnum*ndim1*ndim2;
            double * myqc = qc + threadnum*maxl*ndim2;
            double * mycqc = cqc + threadnum*maxl*maxr;

            double * myqp = myqe;
            if(packed())
                myqp = qp + threadnum*ndim12;

            // read this tensor by Q
            this->ReadByQ(myqp, 1, q);

            if(packed())
            {
                // expand packed matrix
                // ndim1 should equal ndim2
                int index = 0;
                for(int i = 0; i < ndim1; i++)
                    for(int j = 0; j <= i; j++)
                        myqe[i*ndim1+j] = myqe[j*ndim2+i] = myqp[index++];
            }


            // actually do the transformation
            C_DGEMM('T', 'N', lncols, ndim2, ndim1, 1.0, lptr, lncols, myqe, ndim2, 0.0, myqc, ndim2);
            C_DGEMM('N', 'N', lncols, rncols, ndim2, 1.0, myqc, ndim2, rptr, rncols, 0.0, mycqc, rncols);


            // write out
            // Write to memory or disk
            if(qout->packed())
            {
                // can use qc for scratch
                for(int i = 0, index = 0; i < lncols; i++)
                for(int j = 0; j <= i; j++, index++)
                    myqc[index] = mycqc[i*rncols+j];

                qout->WriteByQ_(myqc, 1, q);
            }
            else
                qout->WriteByQ_(mycqc, 1, q);

            #ifdef PANACHE_TIMING
            tim.Stop();
            qout->GenTimer().AddTime(tim);
            #endif
        }
    }

    // done with stuff
    delete [] qe;
    delete [] qc;
    delete [] cqc;

    if(packed())
        delete [] qp;
}



//////////////////////////////
// DiskQTensor
//////////////////////////////
void ThreeIndexTensor::DiskQTensor::OpenFile_(void)
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

void ThreeIndexTensor::DiskQTensor::CloseFile_(void)
{
    if(file_ && file_->is_open())
    {
        file_->close();
        file_.reset();
    }
}

void ThreeIndexTensor::DiskQTensor::Reset_(void)
{
    file_->seekg(0);
    file_->seekp(0);
}

void ThreeIndexTensor::DiskQTensor::Write_(double * data, int nij, int ijstart)
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

void ThreeIndexTensor::DiskQTensor::WriteByQ_(double * data, int nq, int qstart)
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

void ThreeIndexTensor::DiskQTensor::Read_(double * data, int nij, int ijstart)
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

void ThreeIndexTensor::DiskQTensor::ReadByQ_(double * data, int nq, int qstart)
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

void ThreeIndexTensor::DiskQTensor::Clear_(void)
{
    //! \todo Erase file
    CloseFile_();
}

void ThreeIndexTensor::DiskQTensor::Init_(void)
{
    OpenFile_();
}


ThreeIndexTensor::DiskQTensor::DiskQTensor()
{
}



//////////////////////////////
// MemoryQTensor
//////////////////////////////
void ThreeIndexTensor::MemoryQTensor::Reset_(void)
{
    // nothing needed
}


void ThreeIndexTensor::MemoryQTensor::Write_(double * data, int nij, int ijstart)
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


void ThreeIndexTensor::MemoryQTensor::WriteByQ_(double * data, int nq, int qstart)
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

void ThreeIndexTensor::MemoryQTensor::Read_(double * data, int nij, int ijstart)
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

void ThreeIndexTensor::MemoryQTensor::ReadByQ_(double * data, int nq, int qstart)
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
            data[q0*ndim12() + calcindex(m,n)]
                = data_[calcindex(m,n)*naux() + q];
    }
}

void ThreeIndexTensor::MemoryQTensor::Clear_(void)
{
    data_.reset();
}

void ThreeIndexTensor::MemoryQTensor::Init_(void)
{
    if(!data_)
        data_ = std::unique_ptr<double []>(new double[storesize()]);
}

ThreeIndexTensor::MemoryQTensor::MemoryQTensor()
{
}



std::unique_ptr<ThreeIndexTensor::StoredQTensor> 
ThreeIndexTensor::StoredQTensorFactory(int storeflags) const
{
    if(storeflags & QSTORAGE_ONDISK)
        return std::unique_ptr<ThreeIndexTensor::StoredQTensor>(new DiskQTensor());

    #ifdef PANACHE_CYCLOPS
    else if(storeflags & QSTORAGE_CYCLOPS)
        return std::unique_ptr<ThreeIndexTensor::StoredQTensor>(new CyclopsQTensor());
    #endif

    else
        return std::unique_ptr<ThreeIndexTensor::StoredQTensor>(new MemoryQTensor());
}

std::unique_ptr<ThreeIndexTensor::StoredQTensor> 
ThreeIndexTensor::StoredQTensorFactory(int naux, int ndim1, int ndim2,
                                       int storeflags, const std::string & name) const
{
    if(name == "")
        throw RuntimeError("NO NAME SPECIFIED");

    auto ptr = StoredQTensorFactory(storeflags);

    // convert name to filename if on disk
    if(storeflags & QSTORAGE_ONDISK)
    {
        std::string filename(directory_);
        filename.append("/");
        filename.append(name);
        ptr->Init(naux, ndim1, ndim2, storeflags, filename);
    }
    else
        ptr->Init(naux, ndim1, ndim2, storeflags, name);

    return ptr;
}




} // close namespace panache

