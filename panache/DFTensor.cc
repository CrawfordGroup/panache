/*! \file
 * \brief Density fitting tensor generation and manipulation (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/DFTensor.h"
#include "panache/FittingMetric.h"
#include "panache/Output.h"
#include "panache/BasisSet.h"
#include "panache/BasisSetParser.h"
#include "panache/storedqtensor/StoredQTensor.h"
#include "panache/storedqtensor/StoredQTensorFactory.h"

namespace panache {


void DFTensor::PrintHeader_(void) const
{
    output::printf("  ==> LibPANACHE DF Tensor <==\n\n");

    output::printf(" => Primary Basis Set <= \n\n");
    primary_->print_detail();

    output::printf(" => Auxiliary Basis Set <= \n\n");
    auxiliary_->print_detail();
}

void DFTensor::Init_(void)
{
    naux_ = auxiliary_->nbf();
    fittingmetric_ = SharedFittingMetric(new FittingMetric(auxiliary_, nthreads_));
    fittingmetric_->form_eig_inverse();
}

DFTensor::DFTensor(SharedBasisSet primary, SharedBasisSet auxiliary,
                   const std::string & directory,
                   int nthreads) : ThreeIndexTensor(primary, directory, QGEN_DFQSO, nthreads), auxiliary_(auxiliary)
{
    PrintHeader_();
    Init_();
}

DFTensor::DFTensor(SharedBasisSet primary,
                   const std::string & auxpath,
                   const std::string & directory,
                   int nthreads) : ThreeIndexTensor(primary, directory, QGEN_DFQSO, nthreads),
                                   auxiliary_(CreateAuxFromFile_(auxpath, primary->molecule()))
{
    PrintHeader_();
    Init_();
}

UniqueStoredQTensor DFTensor::GenQso(int storeflags) const
{
    // Always gen qso as packed and by q
    auto qso = StoredQTensorFactory(naux_, nso_, nso_, 
                                    storeflags | QSTORAGE_PACKED | QSTORAGE_BYQ, "qso", directory_);

    qso->GenDFQso(fittingmetric_, primary_, auxiliary_, nthreads_);
    return qso;
}


SharedBasisSet DFTensor::CreateAuxFromFile_(const std::string & auxpath, SharedMolecule mol)
{
    // Gaussian input file parser for the auxiliary basis
    SharedBasisSetParser parser(new Gaussian94BasisSetParser);
    return SharedBasisSet(new BasisSet(parser, mol, auxpath));
}

} // close namespace panache
