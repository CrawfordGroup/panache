/*! \file
 * \brief Cholesky tensor generation and manipulation (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/CHTensor.h"
#include "panache/BasisSet.h"
#include "panache/Output.h"
#include "panache/storedqtensor/StoredQTensor.h"
#include "panache/storedqtensor/StoredQTensorFactory.h"

namespace panache {

CHTensor::CHTensor(SharedBasisSet primary, double delta,
                   const std::string & directory,
                   int nthreads) : ThreeIndexTensor(primary, directory, QGEN_CHQSO, nthreads), delta_(delta)
{
    output::printf("  ==> LibPANACHE CH Tensor <==\n\n");
    output::printf("  delta: %f", delta_);

    output::printf(" => Primary Basis Set <= \n\n");
    primary_->print_detail();
}


UniqueStoredQTensor CHTensor::GenQso(int storeflags) const
{
    auto qso = StoredQTensorFactory(storeflags);

    // will be initialized in here
    qso->GenCHQso(primary_, delta_, storeflags, nthreads_);
    return qso;
}

} // close namespace panache
