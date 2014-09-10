#include <fstream>

#include "panache/Exception.h"
#include "panache/DFTensor.h"
#include "panache/MPI.h"

namespace panache
{


//////////////////////////////
// CyclopsQTensor
//////////////////////////////
void DFTensor::CyclopsQTensor::Reset_(void)
{
    // nothing needed
}

void DFTensor::CyclopsQTensor::Write_(double * data, int nij, int ijstart)
{
    size_t nelements = nij*ndim12();

    std::vector<long> indices;
    indices.reserve(nelements);


    IJIterator it(ndim1(), ndim2(), packed());
    it += ijstart;

    if(byq())
    {
        for(int q = 0; q < naux(); q++)
        {
            indices.push_back(q+it.i()*naux()+it.j()*naux()*ndim1());
            ++it;
        }
    }
    else
    {
        for(int q = 0; q < naux(); q++)
        {
            indices.push_back(it.i()+it.j()*ndim1()+q*ndim1()*ndim2());
            ++it;
        }
    }

    tensor_->write(nelements, indices.data(), data);
}

void DFTensor::CyclopsQTensor::WriteByQ_(double * data, int nq, int qstart)
{


}

void DFTensor::CyclopsQTensor::Read_(double * data, int nij, int ijstart)
{
    size_t nelements = nij*ndim12();

    std::vector<long> indices;
    indices.reserve(nelements);


    IJIterator it(ndim1(), ndim2(), packed());
    it += ijstart;

    if(byq())
    {
        for(int q = 0; q < naux(); q++)
        {
            indices.push_back(q+it.i()*naux()+it.j()*naux()*ndim1());
            ++it;
        }
    }
    else
    {
        for(int q = 0; q < naux(); q++)
        {
            indices.push_back(it.i()+it.j()*ndim1()+q*ndim1()*ndim2());
            ++it;
        }
    }

    tensor_->read(nelements, indices.data(), data);
}

void DFTensor::CyclopsQTensor::ReadByQ_(double * data, int nq, int qstart)
{
}

void DFTensor::CyclopsQTensor::Clear_(void)
{
    tensor_.release();
}

void DFTensor::CyclopsQTensor::Init_(void)
{
}


DFTensor::CyclopsQTensor::CyclopsQTensor(int naux, int ndim1, int ndim2, int storeflags, const char *name)
            : StoredQTensor(naux, ndim1, ndim2, storeflags)
{

    int dims[3] = {naux, ndim1, ndim2};

    int syms1[3] = {NS, NS, NS};
    int syms2[3] = {NS, SY, SY};
    int * syms = syms1;

    if(packed())
        syms = syms2;

    tensor_ = std::unique_ptr<CTF_Tensor>(new CTF_Tensor(3, dims, syms, mpi::CTFWorld(), name));
}



void DFTensor::CyclopsQTensor::Transform(const std::vector<TransformMat> & left,
                                         const std::vector<TransformMat> & right,
                                         std::vector<StoredQTensor *> results,
                                         int nthreads)
{
/*
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
    double * qe = new double[ndim1_*ndim2_*nthreads];   // expanded q
    double * qc = new double[maxl*ndim2_*nthreads];     // C(t) Q
    double * cqc = new double[maxl*maxr*nthreads];    // C(t) Q C

    double * qp = qe;
    if(packed())
        qp = new double[ndim12_*nthreads];         // packed q


    #ifdef _OPENMP
        #pragma omp parallel for num_threads(nthreads)
    #endif
    for(int q = 0; q < naux_; q++)
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

            StoredQTensor * qout = results[i];

            double * myqe = qe + threadnum*ndim1_*ndim2_;
            double * myqc = qc + threadnum*maxl*ndim2_;
            double * mycqc = cqc + threadnum*maxl*maxr;

            double * myqp = myqe;
            if(packed())
                myqp = qp + threadnum*ndim12_;

            // read this tensor by Q
            this->ReadByQ(myqp, 1, q);

            if(packed())
            {
                // expand packed matrix
                // ndim1_ should equal ndim2_
                int index = 0;
                for(int i = 0; i < ndim1_; i++)
                    for(int j = 0; j <= i; j++)
                        myqe[i*ndim1_+j] = myqe[j*ndim2_+i] = myqp[index++];
            }


            // actually do the transformation
            C_DGEMM('T', 'N', lncols, ndim2_, ndim1_, 1.0, lptr, lncols, myqe, ndim2_, 0.0, myqc, ndim2_);
            C_DGEMM('N', 'N', lncols, rncols, ndim2_, 1.0, myqc, ndim2_, rptr, rncols, 0.0, mycqc, rncols);


            // write out            
            // Write to memory or disk
            if(qout->packed())
            {
                // can use qc for scratch
                for(int i = 0, index = 0; i < lncols; i++)
                for(int j = 0; j <= i; j++, index++)
                    myqc[index] = mycqc[i*rncols+j];
    
                qout->WriteByQ(myqc, 1, q);
            }
            else
                qout->WriteByQ(mycqc, 1, q);

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
*/
}


} // close namespace panache

