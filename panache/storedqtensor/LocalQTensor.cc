/*! \file
 * \brief Generic, local three-index tensor storage (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <cmath>
#include <fstream>

#include "panache/storedqtensor/LocalQTensor.h"
#include "panache/BasisSet.h"
#include "panache/Lapack.h"
#include "panache/ERI.h"
#include "panache/Flags.h"
#include "panache/Iterator.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace panache
{


//////////////////////////////
// LocalQTensor
//////////////////////////////

LocalQTensor::LocalQTensor(int storeflags, const std::string & name, const std::string & directory) 
              : StoredQTensor(storeflags, name), directory_(directory)
{
    filename_ = directory_;
    filename_.append("/");
    filename_.append(name);
    dimfilename_ = filename_;
    dimfilename_.append(".dim");

    existed_ = FileExists();
}


const std::string & LocalQTensor::filename(void) const { return filename_; }

const std::string & LocalQTensor::directory(void) const { return directory_; }

const std::string & LocalQTensor::dimfilename(void) const { return dimfilename_; }

bool LocalQTensor::FileExists(void) const
{
    std::ifstream ifs(filename_.c_str());
    std::ifstream ifsd(dimfilename_.c_str());
    return (ifs.is_open() && ifsd.is_open());
}


void LocalQTensor::ComputeDiagonal_(std::vector<SharedTwoBodyAOInt> & eris, 
                                                      double * target)
{
    SharedBasisSet basis = eris[0]->basis();

    size_t nthreads = eris.size();


#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
#endif
    for (int M = 0; M < basis->nshell(); M++) 
    {
        int threadnum = 0;

#ifdef _OPENMP
        threadnum = omp_get_thread_num();
#endif
        TwoBodyAOInt * integral = eris[threadnum].get();
        const double* buffer = integral->buffer();

        int nM = basis->shell(M).nfunction();
        int mstart = basis->shell(M).function_index();

        for (int N = 0; N <= M; N++) 
        {
            int nint = integral->compute_shell(M,N,M,N);

            if(nint)
            {
                int nN = basis->shell(N).nfunction();
                int nstart = basis->shell(N).function_index();
    
    
                if(N == M)
                {
                    for (int om = 0; om < nM; om++)
                    for (int on = 0; on <= om; on++)
                        target[((om + mstart) * (om + mstart + 1))/2 + (on + nstart)] =
                            buffer[om * nN * nM * nN + on * nM * nN + om * nN + on];
                }
                else
                {
                    for (int om = 0; om < nM; om++)
                    for (int on = 0; on < nN; on++)
                        target[((om + mstart) * (om + mstart + 1))/2 + (on + nstart)] =
                            buffer[om * nN * nM * nN + on * nM * nN + om * nN + on];
                }
            }
        }
    }
}


void LocalQTensor::ComputeRow_(std::vector<SharedTwoBodyAOInt> & eris, 
                                                 int row, double* target)
{
    SharedBasisSet basis = eris[0]->basis();

    const int nbf = basis->nbf();

    IJIterator ijit(nbf, nbf, true);
    ijit += row; 

    int r = ijit.i();
    int s = ijit.j();
    int R = basis->function_to_shell(r);
    int S = basis->function_to_shell(s);

    int nR = basis->shell(R).nfunction();
    int nS = basis->shell(S).nfunction();
    int rstart = basis->shell(R).function_index();
    int sstart = basis->shell(S).function_index();

    int oR = r - rstart;
    int oS = s - sstart;

    size_t nthreads = eris.size();

    
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
#endif
    for (int M = 0; M < basis->nshell(); M++)
    {
        int threadnum = 0;

#ifdef _OPENMP
        threadnum = omp_get_thread_num();
#endif

        TwoBodyAOInt * integral = eris[threadnum].get();
        const double* buffer = integral->buffer();

        int nM = basis->shell(M).nfunction();
        int mstart = basis->shell(M).function_index();

        for (int N = 0; N <= M; N++)
        {
            int nint = integral->compute_shell(M,N,R,S);

            if(nint)
            {
                int nN = basis->shell(N).nfunction();
                int nstart = basis->shell(N).function_index();
    
                if(N == M)
                {
                    for (int om = 0; om < nM; om++)
                    for (int on = 0; on <= om; on++)
                        target[((om + mstart) * (om + mstart + 1))/2 + (on + nstart)] =
                            buffer[om * nN * nR * nS + on * nR * nS + oR * nS + oS];
                }
                else
                {
                    for (int om = 0; om < nM; om++)
                    for (int on = 0; on < nN; on++)
                        target[((om + mstart) * (om + mstart + 1))/2 + (on + nstart)] =
                            buffer[om * nN * nR * nS + on * nR * nS + oR * nS + oS];
                }
            }
        }
    }
}



void LocalQTensor::GenCHQso_(const SharedBasisSet primary,
                                               double delta,
                                               int nthreads)
{
    // number of threads is passed around implicitly as the size of eris
    std::vector<SharedTwoBodyAOInt> eris;

    for(int i = 0; i < nthreads; i++)
        eris.push_back(GetERI(primary, primary, primary, primary));

    int nQ = 0;
    int n = primary->nbf();
    int n12 = (n*(n+1))/2;

    double * diag = new double[n12];

    // actually important. LibERD interface
    // may not fill every value
    std::fill(diag, diag + n12, 0.0);

    ComputeDiagonal_(eris, diag);

    // Temporary cholesky factor
    std::vector<double*> L;

    // List of selected pivots
    std::vector<int> pivots;
 
    while(nQ < n12)
    {
        int pivot = 0;
        double Dmax = diag[0];
        for(int P = 0; P < n12; P++)
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

        L.push_back(new double[n12]);
        std::fill(L.back(), L.back()+n12, 0.0);

        ComputeRow_(eris, pivot, L[nQ]);

        // [(m|Q) - L_m^P L_Q^P]
        for (int P = 0; P < nQ; P++)
            C_DAXPY(n12,-L[P][pivot],L[P],1,L[nQ],1);

        // 1/L_QQ [(m|Q) - L_m^P L_Q^P]
        C_DSCAL(n12, 1.0 / L_QQ, L[nQ], 1);

        // Zero the upper triangle
        for (size_t P = 0; P < pivots.size(); P++)
            L[nQ][pivots[P]] = 0.0;

        // Set the pivot factor
        L[nQ][pivot] = L_QQ;

        // Update the Schur complement diagonal
        for (int P = 0; P < n12; P++)
            diag[P] -= L[nQ][P] * L[nQ][P];


        // Force truly zero elements to zero
        for (size_t P = 0; P < pivots.size(); P++)
            diag[pivots[P]] = 0.0;

        nQ++;
    }

    delete [] diag;
    pivots.clear();
    eris.clear();

    // copy to memory/disk now that we have the sizes
    StoredQTensor::Init(nQ, n, n);

    for(int i = 0; i < nQ; i++)
    {
        WriteByQ_(L[i], 1, i);
        delete [] L[i];
    }
     
}


void LocalQTensor::Transform_(const std::vector<TransformMat> & left,
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


} // close namespace panache

