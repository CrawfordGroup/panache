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

    // Defaults for fitting metric
    if(optflag_ == 0)
        optflag_ = DFOPT_COULOMB | DFOPT_EIGINV;
}

DFTensor::DFTensor(SharedBasisSet primary, SharedBasisSet auxiliary,
                   const std::string & directory,
                   int optflag, int bsorder, int nthreads) 
           : ThreeIndexTensor(primary, directory, QTYPE_DFQSO, bsorder, nthreads),
             auxiliary_(auxiliary), optflag_(optflag)
{
    PrintHeader_();
    Init_();
}

DFTensor::DFTensor(SharedBasisSet primary,
                   const std::string & auxpath,
                   const std::string & directory,
                   int optflag, int bsorder, int nthreads) 
           : ThreeIndexTensor(primary, directory, QTYPE_DFQSO, bsorder, nthreads),
             auxiliary_(CreateAuxFromFile_(auxpath, primary->molecule())), optflag_(optflag)
{
    PrintHeader_();
    Init_();
}

UniqueStoredQTensor DFTensor::GenQso(int storeflags) const
{
    // Since main options can only be set in the constructor, there is no danger
    // of changing options after construction. Therefore, calculations must
    // be equivalent, except for storage options

    // Always gen qso as packed and by q
    auto qso = StoredQTensorFactory(naux_, nso_, nso_, 
                                    storeflags | QSTORAGE_PACKED | QSTORAGE_BYQ, "qso", directory_);

    // already existed
    if(qso->filled())
        return qso;


    // we only use the fitting metric here. No need to keep it around for longer than
    // is necessary
    // Overhead for shared pointer is negligible here, but consider making it a straight
    // pass by reference in the future
    // NOTE: Keeping it a shared pointer since it may be held by some StoredQTensor
    // derived classes and applied after MO transformation
    SharedFittingMetric fittingmetric(new FittingMetric(auxiliary_, nthreads_));

    if(optflag_ & DFOPT_COULOMB)
        fittingmetric->form_coulomb_fitting_metric();
    else
        throw RuntimeError("Unknown fitting metric type!");
     
    if(optflag_ & DFOPT_EIGINV) 
        fittingmetric->form_eig_inverse();
    else if(optflag_ & DFOPT_CHOINV)
        fittingmetric->form_cholesky_inverse();
    else
        throw RuntimeError("Unknown fitting metric decomposition!");


    qso->GenDFQso(fittingmetric, primary_, auxiliary_, nthreads_);

    fittingmetric.reset(); // done with it?

    return qso;
}


SharedBasisSet DFTensor::CreateAuxFromFile_(const std::string & auxpath, SharedMolecule mol)
{
    // Gaussian input file parser for the auxiliary basis
    SharedBasisSetParser parser(new Gaussian94BasisSetParser);
    return SharedBasisSet(new BasisSet(parser, mol, auxpath));
}

} // close namespace panache
