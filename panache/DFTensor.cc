/*! \file
 * \brief Density fitting tensor generation and manipulation (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include "panache/ERI.h"
#include "panache/DFTensor.h"
#include "panache/FittingMetric.h"
#include "panache/Output.h"
#include "panache/BasisSet.h"
#include "panache/BasisSetParser.h"
#include "panache/storedqtensor/StoredQTensor.h"
#include "panache/storedqtensor/StoredQTensorFactory.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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

void DFTensor::GenQTensors_(int qflags, int storeflags) const
{
    // Since main options can only be set in the constructor, there is no danger
    // of changing options after construction. Therefore, calculations must
    // be equivalent, except for storage options and qflags

    // we only use the fitting metric here. No need to keep it around for longer than
    // is necessary
    UniqueFittingMetric fittingmetric(new FittingMetric(auxiliary_, nthreads_));

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


    // Now we have the metric. Build the three-index tensors 

    // default constructor = zero basis
    SharedBasisSet zero(new BasisSet);

    std::vector<SharedTwoBodyAOInt> eris;
    std::vector<const double *> eribuffers;
    std::vector<double *> buff, buffwork, buff2; // local buffers

    const int maxpershell = auxiliary_->max_function_per_shell();

    for(int i = 0; i < nthreads_; i++)
    {
        eris.push_back(GetERI(auxiliary_, zero, primary_, primary_));
        eribuffers.push_back(eris.back()->buffer());
        buff.push_back(new double[maxpershell*nso2_]);
        buffwork.push_back(new double[maxpershell*nso2_]); // for transformed tensors
        buff2.push_back(new double[maxpershell*nso2_]); // for transformed tensors
    }

    const int nprimshell = primary_->nshell();
    const int nauxshell = auxiliary_->nshell();

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads_)
#endif
    for (int P = 0; P < nauxshell; P++)
    {
        int threadnum = 0;
#ifdef _OPENMP
        threadnum = omp_get_thread_num();
#endif
        double * mybuff = buff[threadnum];
        double * mybuff2 = buff2[threadnum];
        double * mybuffwork = buffwork[threadnum];
        const double * myeribuffer = eribuffers[threadnum];

        int np = auxiliary_->shell(P).nfunction();
        int pstart = auxiliary_->shell(P).function_index();
        int pend = pstart + np;

        for (int M = 0; M < nprimshell; M++)
        {
            int nm = primary_->shell(M).nfunction();
            int mstart = primary_->shell(M).function_index();
            int mend = mstart + nm;

            for (int N = 0; N <= M; N++)
            {
                int nn = primary_->shell(N).nfunction();
                int nstart = primary_->shell(N).function_index();
                int nend = nstart + nn;

                int ncalc = eris[threadnum]->compute_shell(P,0,M,N);

                // keep in mind that we are storing this packed
                // but waiting until after the transformations to pack it
                if(ncalc)
                {
                    for (int p = pstart, p0 = 0; p < pend; p++, p0++)
                    {
                        const int pp0 = p0*nm*nn;
                    
                        for (int m = mstart, m0 = 0; m < mend; m++, m0++)
                        {
                            const int mm = m*nn;
                            const int mm0 = m0*nn;

                            // todo can maybe be a bit more efficient
                            for (int n = nstart, n0 = 0; n < nend; n++, n0++)
                                mybuff[pp0 + mm + n] = mybuff[pp0 + n*nm + m] = myeribuffer[pp0 + mm0 + n0];
                        }
                    }
                }
            }
        }

        // Transform and store these rows if the tensor is calculated
        if(qflags & QGEN_QSO)
        {
            // no transformation needed
            qso_->WriteByQ(mybuff, np, pstart, false);
        }
        if(qflags & QGEN_QMO)
        {
            Transform_(nso_, nso2_, np, mybuff, nmo_, Cmo_.get(), nmo_, Cmo_.get(), mybuff2, mybuffwork);
            qmo_->WriteByQ(mybuff, np, pstart, false);
        }
        if(qflags & QGEN_QOO)
        {
            Transform_(nso_, nso2_, np, mybuff, nocc_, Cmo_occ_.get(), nocc_, Cmo_occ_.get(), mybuff2, mybuffwork);
            qoo_->WriteByQ(mybuff, np, pstart, false);
        }
        if(qflags & QGEN_QOV)
        {
            Transform_(nso_, nso2_, np, mybuff, nocc_, Cmo_occ_.get(), nvir_, Cmo_vir_.get(), mybuff2, mybuffwork);
            qov_->WriteByQ(mybuff, np, pstart, false);
        }
        if(qflags & QGEN_QVV)
        {
            Transform_(nso_, nso2_, np, mybuff, nvir_, Cmo_vir_.get(), nvir_, Cmo_vir_.get(), mybuff2, mybuffwork);
            qvv_->WriteByQ(mybuff, np, pstart, false);
        }
    }

    // done!
    // apply metric


    for(int i = 0; i < nthreads_; i++)
    {
        delete [] buff[i];
        delete [] buffwork[i];
        delete [] buff2[i];
    }

}


SharedBasisSet DFTensor::CreateAuxFromFile_(const std::string & auxpath, SharedMolecule mol)
{
    // Gaussian input file parser for the auxiliary basis
    SharedBasisSetParser parser(new Gaussian94BasisSetParser);
    return SharedBasisSet(new BasisSet(parser, mol, auxpath));
}

} // close namespace panache
