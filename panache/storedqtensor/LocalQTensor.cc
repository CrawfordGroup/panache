/*! \file
 * \brief Generic, local three-index tensor storage (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/storedqtensor/LocalQTensor.h"
#include "panache/BasisSet.h"
#include "panache/FittingMetric.h"
#include "panache/Lapack.h"
#include "panache/ERI.h"
#include "panache/Flags.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace panache
{


//////////////////////////////
// LocalQTensor
//////////////////////////////


LocalQTensor::LocalQTensor()
{
}



void LocalQTensor::GenDFQso_(const std::shared_ptr<FittingMetric> & fit,
                                     const SharedBasisSet primary,
                                     const SharedBasisSet auxiliary,
                                     int nthreads)
{
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
}


void LocalQTensor::ComputeDiagonal_(std::vector<std::shared_ptr<TwoBodyAOInt>> & eris, 
                                                      double * target)
{
    SharedBasisSet basis = eris[0]->basis();

    const int nbf = basis->nbf();

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
    
    
                for (int om = 0; om < nM; om++) {
                    for (int on = 0; on < nN; on++) {
                        target[(om + mstart) * nbf + (on + nstart)] =
                        target[(on + nstart) * nbf + (om + mstart)] =
                            buffer[om * nN * nM * nN + on * nM * nN + om * nN + on];
                    }
                }
            }
        }
    }
}


void LocalQTensor::ComputeRow_(std::vector<std::shared_ptr<TwoBodyAOInt>> & eris, 
                                                 int row, double* target)
{
    SharedBasisSet basis = eris[0]->basis();

    const int nbf = basis->nbf();

    int r = row / nbf;
    int s = row % nbf;
    int R = basis->function_to_shell(r);
    int S = basis->function_to_shell(s);

    int nR = basis->shell(R).nfunction();
    int nS = basis->shell(S).nfunction();
    int rstart = basis->shell(R).function_index();
    int sstart = basis->shell(S).function_index();

    int oR = r - rstart;
    int os = s - sstart;

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
    
                for (int om = 0; om < nM; om++) {
                    for (int on = 0; on < nN; on++) {
                        target[(om + mstart) * nbf + (on + nstart)] =
                        target[(on + nstart) * nbf + (om + mstart)] =
                            buffer[om * nN * nR * nS + on * nR * nS + oR * nS + os];
                    }
                }
            }
        }
    }
}



void LocalQTensor::GenCHQso_(const SharedBasisSet primary,
                                               double delta,
                                               int storeflags,
                                               int nthreads)
{
    std::vector<std::shared_ptr<TwoBodyAOInt>> eris;

    for(int i = 0; i < nthreads; i++)
        eris.push_back(GetERI(primary, primary, primary, primary));

    int nQ = 0;
    int n = primary->nbf();
    int n2 = n*n;

    double * diag = new double[n2];
    std::fill(diag, diag + n2, 0.0);

    ComputeDiagonal_(eris, diag);


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
        std::fill(L.back(), L.back()+n2, 0.0);

        ComputeRow_(eris, pivot, L[nQ]);

        // [(m|Q) - L_m^P L_Q^P]
        for (int P = 0; P < nQ; P++)
            C_DAXPY(n2,-L[P][pivots[nQ]],L[P],1,L[nQ],1);

        // 1/L_QQ [(m|Q) - L_m^P L_Q^P]
        C_DSCAL(n2, 1.0 / L_QQ, L[nQ], 1);

        // Zero the upper triangle
        for (size_t P = 0; P < pivots.size(); P++)
            L[nQ][pivots[P]] = 0.0;

        // Set the pivot factor
        L[nQ][pivot] = L_QQ;

        // Update the Schur complement diagonal
        for (int P = 0; P < n2; P++)
            diag[P] -= L[nQ][P] * L[nQ][P];

        // Force truly zero elements to zero
        for (size_t P = 0; P < pivots.size(); P++)
            diag[pivots[P]] = 0.0;

        nQ++;
    }

    delete [] diag;

    // copy to memory now that we have the sizes
    StoredQTensor::Init(nQ, n, n, storeflags | QSTORAGE_BYQ, "qso");

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

