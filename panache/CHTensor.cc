/*! \file
 * \brief Cholesky tensor generation and manipulation (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/CHTensor.h"
#include "panache/Exception.h"
#include "panache/BasisSet.h"
#include "panache/Output.h"
#include "panache/storedqtensor/StoredQTensor.h"
#include "panache/storedqtensor/StoredQTensorFactory.h"

namespace panache {

CHTensor::CHTensor(SharedBasisSet primary, double delta,
                   const std::string & directory,
                   int bsorder,
                   int nthreads) : ThreeIndexTensor(primary, directory, QTYPE_CHQSO, bsorder, nthreads),
                                   delta_(delta)
{
    output::printf("  ==> LibPANACHE CH Tensor <==\n\n");
    output::printf("  delta: %f", delta_);

    output::printf(" => Primary Basis Set <= \n\n");
    primary_->print_detail();
}


void CHTensor::GenQTensors_(int qflags, int storeflags) const
{
    throw RuntimeError("NYI");
}
/*
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
UniqueStoredQTensor CHTensor::GenQso(int storeflags) const
{
    // Since main options can only be set in the constructor, there is no danger
    // of changing options after construction. Therefore, calculations must
    // be equivalent, except for storage options

    // Can't do full initialization yet. Will be done in GenCHQso (virtual function)
    auto qso = StoredQTensorFactory(storeflags | QSTORAGE_BYQ | QSTORAGE_PACKED, "qso", directory_);

    // already existed
    if(qso->filled())
        return qso;

    // will be initialized in here
    qso->GenCHQso(primary_, delta_, nthreads_);
    return qso;
}
*/

} // close namespace panache
