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


    int C_init(int ncenters,
               C_AtomCenter * atoms,
               int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
               int * aux_nshellspercenter, struct C_ShellInfo * aux_shells)
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


    int C_init2(int ncenters,
                C_AtomCenter * atoms,
                int * primary_nshellspercenter, struct C_ShellInfo * primary_shells,
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




    void C_cleanup(int df_handle)
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

    void C_QAO(int df_handle, double * matout, int matsize)
    {
        // matsize is checked in here
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->Qso(matout, matsize);
    }

    int C_TensorDimensions(int df_handle, int * d1, int * d2, int * d3)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        //DFTensor class takes d1-d3 by reference
        dftensors_[df_handle]->TensorDimensions(*d1, *d2, *d3);
    }


    int C_CalculateERI(int df_handle, double * qso, int qsosize, int shell1, int shell2, int shell3, int shell4, double * outbuffer, int buffersize)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        dftensors_[df_handle]->CalculateERI(qso, qsosize, shell1, shell2, shell3, shell4, outbuffer, buffersize);
    }


    int C_CalculateERIMulti(int df_handle,
                            double * qso, int qsosize,
                            int shell1, int nshell1,
                            int shell2, int nshell2,
                            int shell3, int nshell3,
                            int shell4, int nshell4,
                            double * outbuffer, int buffersize)
    {
        if(dftensors_.count(df_handle) == 0)
            throw RuntimeError("Error - cannot find DFTensor object with that handle!");

        return dftensors_[df_handle]->CalculateERIMulti(qso, qsosize, shell1, nshell1, shell2, nshell2,
                                                        shell3, nshell3, shell4, nshell4, outbuffer, buffersize);
    }

}

