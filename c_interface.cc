#include <map>

#include "c_interface.h"
#include "c_convert.h"
#include "DFTensor.h"

#include "Exception.h"
#include "BasisSetParser.h"

using panache::RuntimeError;

namespace
{

int tensor_index_ = 0;
std::map<int, panache::DFTensor *> dftensors_;

} // close anonymous namespace


extern "C" {


    INTTYPE C_init(INTTYPE ncenters,
               C_AtomCenter * atoms,
               INTTYPE * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
               INTTYPE * aux_nshellspercenter, struct C_ShellInfo * aux_shells)
    {
        // Molecule
        std::shared_ptr<panache::Molecule> molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        auto primaryBasis = panache::BasisSetFromArrays(molecule, ncenters,
                            primary_nshellspercenter, primary_shells);

        auto auxBasis = panache::BasisSetFromArrays(molecule, ncenters,
                        aux_nshellspercenter, aux_shells);

        panache::DFTensor * dft = new panache::DFTensor(primaryBasis, auxBasis);
        dftensors_[tensor_index_] = dft;

        return tensor_index_++;
    }


    INTTYPE C_init2(INTTYPE ncenters,
                    C_AtomCenter * atoms,
                    INTTYPE * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
                    const char * auxfilename)
    {
        // Molecule
        std::shared_ptr<panache::Molecule> molecule = panache::MoleculeFromArrays(ncenters, atoms);


        // Construct the basis set info
        auto primaryBasis = panache::BasisSetFromArrays(molecule, ncenters,
                            primary_nshellspercenter, primary_shells);

        // Gaussian input file parser for the auxiliary basis
        std::shared_ptr<panache::Gaussian94BasisSetParser> parser(new panache::Gaussian94BasisSetParser);
        auto auxBasis = panache::BasisSet::construct(parser, molecule, auxfilename);

        panache::DFTensor * dft = new panache::DFTensor(primaryBasis, auxBasis);
        dftensors_[tensor_index_] = dft;

        return tensor_index_++;
    }




    void C_cleanup(INTTYPE df_handle)
    {
        if(dftensors_.count(df_handle) > 0)
        {
            delete dftensors_[df_handle];
            dftensors_.erase(df_handle);
        }
        else
            throw RuntimeError("Error - cannot erase DFTensor object with that handle!");

    }



    void C_cleanup_all(void)
    {
        for(auto it : dftensors_)
            delete it.second;

        dftensors_.clear();
    }

    void C_QAO(INTTYPE df_handle, double * matout, INTTYPE matsize)
    {
        // matsize is checked in here
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->Qso(matout, matsize);
    }

    INTTYPE C_TensorDimensions(INTTYPE df_handle, INTTYPE * d1, INTTYPE * d2, INTTYPE * d3)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        //DFTensor class takes d1-d3 by reference
        // use regular ints first
        int t1, t2, t3;
        dftensors_[df_handle]->TensorDimensions(t1, t2, t3);

        *d1 = t1;
        *d2 = t2;
        *d3 = t3;
    }


    INTTYPE C_CalculateERI(INTTYPE df_handle, double * qso, INTTYPE qsosize, INTTYPE shell1, INTTYPE shell2, INTTYPE shell3, INTTYPE shell4, double * outbuffer, INTTYPE buffersize)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->CalculateERI(qso, qsosize, shell1, shell2, shell3, shell4, outbuffer, buffersize);
    }


    INTTYPE C_CalculateERIMulti(INTTYPE df_handle,
                            double * qso, INTTYPE qsosize,
                            INTTYPE shell1, INTTYPE nshell1,
                            INTTYPE shell2, INTTYPE nshell2,
                            INTTYPE shell3, INTTYPE nshell3,
                            INTTYPE shell4, INTTYPE nshell4,
                            double * outbuffer, INTTYPE buffersize)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        return dftensors_[df_handle]->CalculateERIMulti(qso, qsosize, shell1, nshell1, shell2, nshell2,
                                                        shell3, nshell3, shell4, nshell4, outbuffer, buffersize);
    }

    void C_ReorderQ_GAMESS(INTTYPE df_handle, double * qso, INTTYPE qsosize)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->ReorderQ_GAMESS(qso, qsosize);
    }
}

