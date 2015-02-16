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

/*
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
*/


void LocalQTensor::Transform_(int nleft, double * left,
                              int nright, double * right,
                              int nthreads)
{
    int naux = StoredQTensor::naux();
    int ndim1 = StoredQTensor::ndim1();
    int ndim2 = StoredQTensor::ndim2();
    int ndim12 = StoredQTensor::ndim12();

    // temporary space
    double * qe = new double[ndim1*ndim2*nthreads];   // expanded q
    double * qc = new double[nleft*ndim2*nthreads];    // C(t) Q
    double * cqc = new double[nleft*nright*nthreads];

    double * qp = cqc;
    if(packed())
        qp = new double[ndim12*nthreads];

    #ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads)
    #endif
    for(int q = 0; q < naux; q++)
    {
        int threadnum = 0;

        #ifdef _OPENMP
            threadnum = omp_get_thread_num();
        #endif

        double * myqe = qe + threadnum*ndim1*ndim2;
        double * myqc = qc + threadnum*nleft*ndim2;
        double * mycqc = cqc + threadnum*nleft*nright;
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
        C_DGEMM('T', 'N', nleft, ndim2, ndim1, 1.0, left, nleft, myqe, ndim2, 0.0, myqc, ndim2);
        C_DGEMM('N', 'N', nright, nright, ndim2, 1.0, myqc, ndim2, right, nright, 0.0, mycqc, nright);


        // write out
        // Write to memory or disk
        this->WriteByQ_(mycqc, 1, q, false);
    }

    // done with stuff
    delete [] qe;
    delete [] qc;
    delete [] cqc;

    if(packed())
        delete [] qp;
}

void LocalQTensor::ContractQ_(int n, double * mat, int nthreads)
{
    //! \todo could be done in batches
    const int indim12 = ndim12();
    const int inaux = naux();

    double * work = new double[inaux * nthreads]; 
    double * work2 = new double[inaux * nthreads]; 

    #ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads)
    #endif
    for(int ij = 0; ij < indim12; ij++)
    {
        #ifdef _OPENMP
            const int threadnum = omp_get_thread_num();
        #else
            const int threadnum = 0;
        #endif

        double * mywork = work + threadnum*inaux;
        double * mywork2 = work2 + threadnum*inaux;

        // get the data
        Read(mywork, 1, ij); 

        // multiply (n x naux  with naux x 1)
        C_DGEMM('N', 'N', n, 1, inaux, 1.0, mat, inaux, mywork, 1, 0.0, mywork2, 1);

        // store the data
        Write(mywork2, 1, ij); 
    }

    delete [] work;
    delete [] work2;
}

} // close namespace panache

